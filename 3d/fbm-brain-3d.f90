      PROGRAM fbm3d_brain_6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   time-discrete reflected fractional Brownian motion
!   in 3D shape read as a series of cross sections from files
!
!   15 Jun 2019      fbm_xy_lr7       variable number of holes, reads ymax, ymin   
!   17 Apr 2022      fbm3d_brain_1    first 3d version 
!   10 May 2022      fbm3d_brain_2    add management of holes    
!   25 May 2022      fbm3d_brain_3    integer GRIDMAX for array dimensions
!   03 Jun 2022      fbm3d_brain_4    several allowed areas per section
!   06 Jun 2022      fbm3d_brain_5    add exponential tempering
!   10 Jun 2022      fbm3d_brain_6    specify starting point (cell bodies)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preprocessor directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!define PARALLEL
#define VERSION 'fbm3d_brain_6'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      integer(i4b), parameter     :: M=23,NT=2**M               ! number of time steps (in which mu const) 
      integer(i4b), parameter     :: NCONF=12000                  ! number of walkers
      real(r8b), parameter        :: GAMMA = 0.4D0              ! FBM correlation exponent 
      logical, parameter          :: FRACTIONAL=.true.           ! true: FBM, false: BM

      real(r8b),parameter         :: X0=100.D0, Y0=100.D0, Z0=20.D0     ! starting point
	  real(r8b),parameter         :: X0MIN=-50.D0, X0MAX=50.D0          ! starting point range 
	  real(r8b),parameter         :: Y0MIN=575.D0, Y0MAX=750.D0          ! starting point range 
	  real(r8b),parameter         :: Z0MIN=540.D0, Z0MAX=600.D0          ! starting point range 
	  
	  
      logical, parameter          :: RANDOMSTART=.true.          ! if true start from random position

      real(r8b), parameter        :: STEPSIG=0.4D0                ! sigma of individual step 
      character(3), parameter     :: STEPDIS='GAU'                 ! Gaussian = GAU, binary = BIN, box = BOX                   
	  integer(i4b), parameter     :: NTTEMP=2000                       ! tempering time
	  logical, parameter          :: TEMPERING=.false.              ! switches on tempering of the FBM 
 
	  
	  integer(i4b), parameter     :: SECSTART=1,SECEND=50          ! first and last section to be read
	  real(r8b), parameter        :: SECPITCH = 12.D0              ! z-distance between sections
	  
	  integer(i4b), parameter     :: MAXSHAPES = 2                ! maximum number of outer shapes in section 	  
	  integer(i4b), parameter     :: MAXHOLES = 5                  ! maximum number of holes in sections  

      real(r8b), parameter        :: GRIDXMIN=-391.D0, GRIDXMAX=391.D0   ! dimensions of grid covering brain
      real(r8b), parameter        :: GRIDYMIN=418.D0, GRIDYMAX=920.D0
      real(r8b), parameter        :: GRIDZMIN=10.D0, GRIDZMAX=602.D0     
	  real(r8b), parameter        :: CELLSIZE=2.D0                      ! size of grid cell  
       
      
      logical, parameter          :: WRITETRAJEC = .false.        ! write individual trajectories
      logical, parameter          :: WRITEDISTRIB = .true.       ! write density distribution
      
      integer(i4b), parameter     :: IRINIT=1                    ! random number seed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      integer(i4b), parameter     :: YYMIN=NINT(GRIDYMIN), YYMAX=NINT(GRIDYMAX)
      integer(i4b), parameter     :: NGX=NINT((GRIDXMAX-GRIDXMIN)/CELLSIZE)          ! number of grid cells in x direction
      integer(i4b), parameter     :: NGY=NINT((GRIDYMAX-GRIDYMIN)/CELLSIZE)          ! number of grid cells in y direction
      integer(i4b), parameter     :: NGZ=NINT((GRIDZMAX-GRIDZMIN)/CELLSIZE)          ! number of grid cells in z direction
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(r8b)              :: xx(0:NT),yy(0:NT),zz(0:NT)          ! walker coordinates
      real(r8b)              :: xix(1:NT),xiy(1:NT),xiz(1:NT)       ! increments
      real(r8b)              :: rad2                                ! distance of walker from origin
      
	  integer(i4b)           :: shapelist(MAXSHAPES,SECSTART:SECEND)                                    ! list of which shapes exist (1) or not (0) 
	  integer(i4b)           :: ymin(MAXSHAPES,SECSTART:SECEND), ymax(MAXSHAPES,SECSTART:SECEND)            !  ymin, ymax for each section outline
      integer(i4b)           :: xleft(MAXSHAPES,SECSTART:SECEND,YYMIN:YYMAX),xright(MAXSHAPES,SECSTART:SECEND,YYMIN:YYMAX)     !  left and right boundaries for each section 
	  
	  integer(i4b)           :: holelist(MAXHOLES,SECSTART:SECEND)                                    ! list of which holes exist (1) or not (0) 
	  integer(i4b)           :: yhmin(MAXHOLES,SECSTART:SECEND), yhmax(MAXHOLES, SECSTART:SECEND)     !  ymin, ymax for each hole
      integer(i4b)           :: xhleft(MAXHOLES,SECSTART:SECEND,YYMIN:YYMAX), xhright(MAXHOLES,SECSTART:SECEND,YYMIN:YYMAX)     !  left and right hole boundaries 
	        
      logical                :: outside
      
      real(r8b)              :: confxx(1:NT)                     ! average of xx after each time step  
      real(r8b)              :: conf2xx(1:NT)                    
      real(r8b)              :: sumxx(1:NT),sum2xx(1:NT)         ! sums over machines in MPI version
      real(r8b)              :: auxxx(1:NT),aux2xx(1:NT) 

      integer(i4b)           :: iconf, isec, it, iy, ishape, ihole             ! configuration, section, time, shape, and hole counters   
      integer(i4b)           :: totconf                         ! actual number of confs
      
      real(r8b)              :: confdens(NGX,NGY,NGZ)                ! density array
      real(r8b)              :: normdens, lognormdens             ! normalized density 
      real(r8b)              :: cellx,celly,cellz                      ! center coordinates of each cell
      integer(i4b)           :: igx,igy,igz                          ! grid cell counter
      real(r8b)              :: sumdens(NGX,NGY,NGZ)                ! sums over machines in MPI version
      real(r8b)              :: auxdens(NGX,NGY,NGZ)                ! 
      
      external               :: kissinit 
      real(r8b),external     :: rkiss05,gkiss05, erfcc 

      character(13)          :: shapefile =   "shapXX_00.dat"
      character(13)          :: holefile =    "holeXX_00.dat"       ! hole number XX
	  character(13)          :: shapelistfile ="shaplist.dat" 
	  character(12)          :: holelistfile ="holelist.dat" 
     
      character(13)          :: trafile = 'tra000000.dat'
      character(11)          :: disfile = 'distrib.dat'

      integer(i4b)           :: dummy
      character(6)           :: cdummy 
            
! Now the MPI stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      include 'mpif.h'
      integer(i4b)              :: ierr
      integer(i4b)              :: id,myid                  ! process index
      integer(i4b)              :: numprocs              ! total number of processes
      integer(i4b)              :: status(MPI_STATUS_SIZE)      
#endif             
  

! Start of main program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set up MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
      call MPI_INIT(ierr)
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      totconf=(NCONF/numprocs)*numprocs
      
      if (myid==0) then
         print *,'Program ',VERSION,' runing on', numprocs, ' processes'
         print *,'--------------------------------------------------'
      endif ! of if (myid==0)
#else
      totconf=NCONF
      print *,'Program ',VERSION,' runing on single processor'
      print *,'--------------------------------------------------'
#endif 

! Read geometry information - shape oulines
      open(8,file=shapelistfile,status='old')
      rewind(8)
      read(8,*) cdummy
      read(8,*) cdummy
	  do isec = SECSTART,SECEND
		 read(8,*) cdummy, shapelist(:,isec)
!         print *, isec, shapelist(:,isec)
	  enddo
	  close(8)

      do isec = SECSTART,SECEND
	  do ishape = 1, MAXSHAPES
         write(shapefile(8:9),'(I2.2)') isec    
         write(shapefile(5:6),'(I2.2)') ishape
		 if (shapelist(ishape,isec).eq.1) then		 
		   open(8,file=shapefile,status='old')
           rewind(8)
           read(8,*) cdummy, ymin(ishape,isec)
           read(8,*) cdummy, ymax(ishape,isec)
!           print *, ymin(ishape,isec),ymax(ishape,isec)
           do iy=ymin(ishape,isec),ymax(ishape,isec)
             read (8,*) dummy, xleft(ishape,isec,iy),xright(ishape,isec,iy)
           enddo
           close(8) 
		 endif
      enddo		 
      enddo

! Read geometry information - holes
      open(8,file=holelistfile,status='old')
      rewind(8)
      read(8,*) cdummy
      read(8,*) cdummy
	  do isec = SECSTART,SECEND
		 read(8,*) cdummy, holelist(:,isec)
!         print *, isec, holelist(:,isec)
	  enddo
	  close(8)

      do isec = SECSTART,SECEND
	  do ihole = 1, MAXHOLES
         write(holefile(8:9),'(I2.2)') isec    
		 write(holefile(5:6),'(I2.2)') ihole
         if (holelist(ihole,isec).eq.1) then		 
           open(8,file=holefile,status='old')
           rewind(8)
           read(8,*) cdummy, yhmin(ihole,isec)
           read(8,*) cdummy, yhmax(ihole,isec)
!           print *, isec, ihole, yhmin(ihole,isec),yhmax(ihole,isec)
           do iy=yhmin(ihole,isec),yhmax(ihole,isec)
              read (8,*) dummy, xhleft(ihole,isec,iy),xhright(ihole,isec,iy)
           enddo
           close(8) 
		 endif
      enddo
	  enddo
	  
      
      confxx(:)=0.D0 
      conf2xx(:)=0.D0 

      confdens(:,:,:0)=0.D0
      
! Loop over disorder configurations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      disorder_loop: do iconf=myid+1,totconf,numprocs
!      if (myid==0) print *, 'dis. conf.', iconf
#else
      disorder_loop: do iconf=1,totconf
         print *, 'dis. conf.', iconf
#endif 
        call gkissinit(IRINIT+iconf-1)

        if (FRACTIONAL) then
          call corvec(xix,NT,M,GAMMA,TEMPERING,NTTEMP)                         ! create x-increments (correlated Gaussian random numbers)
          call corvec(xiy,NT,M,GAMMA,TEMPERING,NTTEMP)                         ! create y-increments (correlated Gaussian random numbers)
          call corvec(xiz,NT,M,GAMMA,TEMPERING,NTTEMP)                         ! create y-increments (correlated Gaussian random numbers)
        else
          do it=1,NT  
            xix(it)=gkiss05()                                 ! create uncorrelated random step
            xiy(it)=gkiss05()                                 ! create uncorrelated random step
            xiz(it)=gkiss05()                                 ! create uncorrelated random step
          enddo 
        endif  
                
        if (STEPDIS.eq.'BOX') then
          xix(:) = 1 - (0.5D0*erfcc(xix(:)/sqrt(2.0D0)))                       ! map onto unit inteval
          xix(:)=xix(:)-0.5D0                                                  ! center around 0
          xiy(:) = 1 - (0.5D0*erfcc(xiy(:)/sqrt(2.0D0)))                       ! map onto unit inteval
          xiy(:)=xiy(:)-0.5D0                                                  ! center around 0
          xiz(:) = 1 - (0.5D0*erfcc(xiz(:)/sqrt(2.0D0)))                       ! map onto unit inteval
          xiz(:)=xiz(:)-0.5D0                                                  ! center around 0

        endif
        
        if (STEPDIS.eq.'BIN') then 
          xix(:)=sign(1.D0,xix(:))
          xiy(:)=sign(1.D0,xiy(:))
          xiz(:)=sign(1.D0,xiz(:))
        endif  
        
        xix(:)=xix(:)*STEPSIG
        xiy(:)=xiy(:)*STEPSIG
        xiz(:)=xiz(:)*STEPSIG
		
        
! Time loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           call initialize_point 
           do it=1, NT/2
              xx(it)=xx(it-1)+xix(it)
              yy(it)=yy(it-1)+xiy(it)
              zz(it)=zz(it-1)+xiz(it)
      
              if (pointnotvalid(xx(it),yy(it),zz(it))) then
                 xx(it)=xx(it-1)
                 yy(it)=yy(it-1)
                 zz(it)=zz(it-1)
              end if
              
              igx= 1+int(xx(it)-GRIDXMIN)/CELLSIZE
              igy= 1+int(yy(it)-GRIDYMIN)/CELLSIZE
              igz= 1+int(zz(it)-GRIDZMIN)/CELLSIZE
			  if ((igz.ge.1).and.(igz.le.NGZ)) then
                confdens(igx,igy,igz)=confdens(igx,igy,igz)+1.D0
			  endif	
           end do  

        if (WRITETRAJEC) then   
          write(trafile(4:9),'(I6.6)') IRINIT+iconf-1       
          open(2,file=trafile,status='replace')
          write(2,*) 'Program ', VERSION
          write(2,*) 'trajectory of reflected FBM on 2d shape'
          write(2,*) 'step distribution ',STEPDIS
          write(2,*) 'step size ',STEPSIG
          write(2,*) 'NT= ', NT
          write(2,*) 'FRACTIONAL=', FRACTIONAL
          write(2,*) 'GAMMMA=', GAMMA
		  write(2,*) 'TEMPERING=', TEMPERING
		  write(2,*) 'NTTEMP=', NTTEMP
          write(2,*) 'IRINIT=',IRINIT
          write(2,*) 'RNG seed=', IRINIT+iconf-1
          write (2,*)'=================================='
          write(2,*) '   time         x         y        z'
          do it=0, NT/2
            Write(2,'(1X,I7,6(2X,E13.6))')  it, xx(it), yy(it), zz(it)
          enddo 
          close(2) 
        endif  

      end do disorder_loop      ! of do incof=1,NCONF
      
! Now collect and communicate data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
     if (myid.ne.0) then                                                  ! Send data
         call MPI_SEND(confdens,NGX*NGY*NGZ,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
      else
         sumdens(:,:,:)=confdens(:,:,:)
         do id=1,numprocs-1                                                   ! Receive data
            call MPI_RECV(auxdens,NGX*NGY*NGZ,MPI_DOUBLE_PRECISION,id,1,MPI_COMM_WORLD,status,ierr)
            sumdens(:,:,:)=sumdens(:,:,:)+auxdens(:,:,:) 
         enddo
      endif        
#else          
      sumdens(:,:,:)=confdens(:,:,:)
#endif       
      
#ifdef PARALLEL
      if (myid==0) then
#endif
      if (WRITEDISTRIB) then
          open(2,file=disfile,status='replace')
#ifdef PARALLEL          
          write(2,*) 'Program ', VERSION, ' running on ', numprocs, ' CPUs'
#else          
          write(2,*) 'Program ', VERSION, ' running serially'
#endif
          write(2,*) 'density distribution of reflected FBM in 3d shape'
          write(2,*) 'step distribution ',STEPDIS
          write(2,*) 'step size ',STEPSIG
          write(2,*) 'NT= ', NT
          write(2,*) 'FRACTIONAL=', FRACTIONAL
          write(2,*) 'GAMMMA=', GAMMA
		  write(2,*) 'TEMPERING=', TEMPERING
		  write(2,*) 'NTTEMP=', NTTEMP
          write(2,*) 'Number of walkers', totconf
          write(2,*) 'IRINIT=',IRINIT
          write (2,*)'=================================='
          write(2,*) '   igx  igy  igz    x      y   z   distrib'
          do igx=1,NGX
          do igy=1,NGY
!		  do igz=1,NGZ,4 
		  do igz=1,NGZ,1 
            cellx= GRIDXMIN + (1.D0*igx -0.5D0)*CELLSIZE
            celly= GRIDYMIN + (1.D0*igy -0.5D0)*CELLSIZE
			cellz= GRIDZMIN + (1.D0*igz -0.5D0)*CELLSIZE
            normdens=sumdens(igx,igy,igz)/totconf/(NT/2)/CELLSIZE**3
            if (normdens.gt.0.D0) then
               lognormdens=log(normdens)
            else
               lognormdens=-100
            endif   
            Write(2,'(3(1X,I4),3(1X,F8.1),6(1X,E12.5))')  igx, igy, igz, cellx, celly, cellz, normdens, lognormdens
          enddo 
          enddo
		  enddo
          close(2)             
      endif
#ifdef PARALLEL
      endif ! of if (myid==0)
#endif        

      
#ifdef PARALLEL
      call MPI_FINALIZE(ierr)
#endif 
      stop      
      
! Now the internal subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
            
      function pointnotvalid(xn,yn,zn)
      real(r8b) :: xn,yn,zn
	  integer(i4b) :: secbef,secaft                   ! sections before and after point in z-direction
	  integer(i4b) :: ybef,yaft                        ! indices corresponding to yn in the before and after sections 
	  real(r8b) :: yfactor,zfactor                      ! scale factors for interpolation
	  integer(i4b) :: ymini, ymaxi                       ! interpolated ymin and ymax 
	  real(r8b)    :: xlefti,xrighti                    ! interpolated left and right boundaries
      logical   :: inshape(MAXSHAPES)                   ! is point in shape 
      logical   :: pointnotvalid
	   
	  
! Test z-direction      
	  if ( (zn.ge.SECEND*SECPITCH).or. (zn.le.SECSTART*SECPITCH)  ) then                  ! outside in z-direction
	     pointnotvalid=.true. 
		 return
      endif
         	  
! Find ajecent sections	  
      secbef = int(zn/SECPITCH)
	  secaft = secbef+1
	  zfactor= (zn-secbef*SECPITCH)/SECPITCH
	  if ( (secbef.lt.SECSTART).or.(secaft.gt.SECEND) ) then
	     print *,'invalid section selected'
         stop
      endif		 

! Test whether point is in one of the shape outlines
      pointnotvalid=.true.
      do ishape = 1, MAXSHAPES	    
	    if ( (shapelist(ishape,secbef).eq.1).and.(shapelist(ishape,secaft).eq.1) ) then
! Interpolate ymin and ymax of shape outline
          ymini = nint( zfactor*ymin(ishape,secaft) + (1.D0-zfactor)*ymin(ishape,secbef) )
	      ymaxi = nint( zfactor*ymax(ishape,secaft) + (1.D0-zfactor)*ymax(ishape,secbef) )
          if ( (nint(yn).ge.ymini).and.(nint(yn).le.ymaxi) ) then
! Interpolate left and right boundaries of outline
            yfactor = dble(nint(yn)-ymini)/dble(ymaxi-ymini)
	        ybef=nint( ymin(ishape,secbef)+yfactor*(ymax(ishape,secbef) - ymin(ishape,secbef)) )
            yaft=nint( ymin(ishape,secaft)+yfactor*(ymax(ishape,secaft) - ymin(ishape,secaft)) ) 
	        xlefti =  zfactor*xleft(ishape,secaft,yaft) + (1.D0-zfactor)*xleft(ishape,secbef,ybef) 
	        xrighti = zfactor*xright(ishape,secaft,yaft) + (1.D0-zfactor)*xright(ishape,secbef,ybef) 
	        if ( (abs(xn).ge.xlefti).and.(abs(xn).le.xrighti) ) then                           ! inside of shape
	          pointnotvalid=.false.
		    endif
		  endif
	    endif	 
	  enddo
      if (pointnotvalid) return

	  
! Test whether point is in a hole
      do ihole = 1, MAXHOLES	  
	    if ( (holelist(ihole,secbef).eq.1).and. (holelist(ihole,secaft).eq.1) ) then     ! hole exists in sections before and after point
! Interpolate ymin and ymax of outline
          ymini = nint( zfactor*yhmin(ihole,secaft) + (1.D0-zfactor)*yhmin(ihole,secbef) )
	      ymaxi = nint( zfactor*yhmax(ihole,secaft) + (1.D0-zfactor)*yhmax(ihole,secbef) )
          if ( (nint(yn).gt.ymini).and.(nint(yn).lt.ymaxi) ) then                         ! point in y-range of hole
! Interpolate left and right boundaries of outline
            yfactor = dble(nint(yn)-ymini)/dble(ymaxi-ymini)
	        ybef=nint( yhmin(ihole,secbef)+yfactor*(yhmax(ihole,secbef) - yhmin(ihole,secbef)) )
            yaft=nint( yhmin(ihole,secaft)+yfactor*(yhmax(ihole,secaft) - yhmin(ihole,secaft)) ) 
	        xlefti =  zfactor*xhleft(ihole,secaft,yaft) + (1.D0-zfactor)*xhleft(ihole,secbef,ybef) 
	        xrighti = zfactor*xhright(ihole,secaft,yaft) + (1.D0-zfactor)*xhright(ihole,secbef,ybef) 
            if ( (abs(xn).gt.xlefti).and.(abs(xn).lt.xrighti) ) then
	          pointnotvalid=.true.
	          return                                                                     ! point in hole
            endif 			  
          endif    
	    endif
	  enddo
	   
	  return
      end function pointnotvalid


      subroutine initialize_point
      if (RANDOMSTART) then            
         xx(0)=0.D0
         yy(0)=0.D0  
		 zz(0)=0.D0
         do while (pointnotvalid(xx(0),yy(0),zz(0)))
           xx(0)=X0MIN+rkiss05()*(X0MAX-X0MIN) 
           yy(0)=Y0MIN+rkiss05()*(Y0MAX-Y0MIN) 
           zz(0)=Z0MIN+rkiss05()*(Z0MAX-Z0MIN) 
         enddo
      else
         xx(0)=X0
         yy(0)=Y0  
		 yy(0)=Z0
      endif   

      end subroutine initialize_point
            
      END PROGRAM fbm3d_brain_6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine CORVEC(xr,Ns,M,gam,tempering,nstemp)
!
! generates a 1d array of Ns Gaussian random numbers xr
! correlated according to a (translationally invariant)
! user-supplied correlation function corfunc(is,Ns)
! 
! uses Fourier filtering method
!
! history
!      v0.9         Dec  7, 2013:        first version, uses Tao Pang FFT
!      v0.91        Oct 11, 2017:        uses much faster FFT by Ooura (in fftsg.f)        
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE corvec(xr,Ns,M,gam,tempering,nstemp)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

      integer(i4b)           :: Ns              ! number of sites, must be power of 2 
      integer(i4b)           :: M               ! Ns=2^M

      real(r8b)              :: xr(0:Ns-1)      ! random number array  
      real(r8b)              :: cr(0:Ns-1)      ! correlation function 
      integer(i4b)           :: is
      
      integer(i4b)           :: ip(0:2+sqrt(1.*Ns))   ! workspace for FFT code
      real(r8b)              :: w(0:Ns/2-1)           ! workspace for FFT code 

      real(r8b)              :: gam                   ! FBM exponent, pass through to correlation function   
	  logical                :: tempering             ! pass through to  correlation function      
      integer(i4b)           :: nstemp                ! tempering time, pass through to  correlation function      
      
      real(r8b), external    :: gkiss05,erfcc
      external               :: rdft                   ! from Ooura's FFT package
      real(r8b),external     :: fbmcorfunc 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if (Ns.ne.2**M) STOP 'Size indices do not match'
! Calculate correlation function 
      do is=0,Ns-1 
         cr(is)= fbmcorfunc(is,Ns,gam,tempering,nstemp) 
      enddo
! Real FFT of correlation function
       ip(0)=0
       call rdft(Ns, 1, cr, ip, w)  
       
! Create array of independent Gaussian random numbers
      do is=0,Ns-1
         xr(is)= gkiss05()
      enddo
! Real FFT of input random numbers      
       call rdft(Ns, 1, xr, ip, w)  
! Filter the Fourier components of the random numbers
! as real-space correlations are symmmetric, FT of c is real
      do is=1,Ns/2-1
          xr(2*is)=xr(2*is)*sqrt(abs(cr(2*is)))*2.D0/Ns
          xr(2*is+1)=xr(2*is+1)*sqrt(abs(cr(2*is)))*2.D0/Ns
      enddo
      xr(0)=xr(0)*sqrt(abs(cr(0)))*2.D0/Ns
      xr(1)=xr(1)*sqrt(abs(cr(1)))*2.D0/Ns 
      
! FFT of filtrered random numbers (back to real space)
       call rdft(Ns, -1, xr, ip, w)  
       
! Transform from Gaussian distribution to flat distribution on (0,1)      
!      do is = 0,Ns-1
!        xr(is) = 1 -  0.5D0*erfcc(xr(is)/sqrt(2.0D0)) 
!      end do
      
      return
      END SUBROUTINE corvec 
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FBM correlation function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION fbmcorfunc(i,N,gam,tempering,nstemp)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      real(r8b)              :: corfunc,fbmcorfunc
      real(r8b)              :: gam
      integer(i4b)           :: i,N
      integer(i4b)           :: dist
	  logical                :: tempering
      integer(i4b)           :: nstemp                       ! tempering time 	  

!      print *, 'gamma=',gam
      dist=min(i,N-i)
      if (dist.eq.0) then 
         corfunc=1.D0
      elseif (dist.lt.1000) then   
         corfunc = 0.5D0*( dble(dist+1)**(2.D0-gam) - 2.D0*(dble(dist)**(2.D0-gam)) + dble(dist-1)**(2.D0-gam) )
      else 
         corfunc = (2.D0-gam)*(1.D0-gam)*dble(dist)**(-gam) 
         corfunc = corfunc + (2.D0-gam)*(1.D0-gam)*(-gam)*(-1.D0-gam)*dble(dist)**(-gam-2.D0)/12.D0          
         corfunc = 0.5D0*corfunc 
      endif   

      if (tempering) then
	     fbmcorfunc=corfunc  * exp( -dble(dist)/dble(nstemp) )
      else
         fbmcorfunc=corfunc
      endif 
      
      return
      END FUNCTION fbmcorfunc 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gaussian random number generator gkiss05
!
! generates normally distributed independent random numbers
! (zero mean, variance 1) using Box-Muller method in polar form
!
! uniform random numbers provided by Marsaglia's kiss (2005 version)
!
! before using the RNG, call gkissinit(seed) to initialize
! the generator. Seed should be a positive integer.
!
!
! History:
!      v0.9     Dec  6, 2013:   first version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION gkiss05()
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      real(r8b)             :: gkiss05
      real(r8b), external   :: rkiss05

      real(r8b)             :: v1,v2,s,fac
      integer(i4b)          :: iset               ! switches between members of the Box-Muller pair
      real(r8b)             :: gset
      common /gausscom/gset,iset

      if (iset.ne.1) then
        do
          v1 = 2.D0 * rkiss05() - 1.D0
          v2 = 2.D0 * rkiss05() - 1.D0
          s = v1 * v1 + v2 * v2
          if ((s<1.D0) .and. (s>0.D0)) exit
        enddo
! Box-Muller transformation creates pairs of random numbers
        fac = sqrt(-2.D0 * log(s) / s)
        gset = v1 * fac
        iset = 1
        gkiss05 = v2 * fac
      else
        iset = 0
        gkiss05 = gset
      end if
      return
      END FUNCTION gkiss05


      SUBROUTINE gkissinit(iinit)
      implicit none
      integer,parameter     :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter     :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      integer(i4b)          :: iinit,iset
      real(r8b)             :: gset
      common /gausscom/gset,iset

      iset=0                         ! resets the switch between the members of the Box-Muller pair
      call kissinit(iinit)           ! initializes the rkiss05 RNG
      end subroutine gkissinit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Random number generator KISS05 after a suggestion by George Marsaglia
! in "Random numbers for C: The END?" posted on sci.crypt.random-numbers
! in 1999
!
! version as in "double precision RNGs" in  sci.math.num-analysis
! http://sci.tech-archive.net/Archive/sci.math.num-analysis/2005-11/msg00352.html
!
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
! Overall period > 2^123
!
!
! A call to rkiss05() gives one random real in the interval [0,1),
! i.e., 0 <= rkiss05 < 1
!
! Before using rkiss05 call kissinit(seed) to initialize
! the generator by random integers produced by Park/Millers
! minimal standard LCG.
! Seed should be any positive integer.
!
! FORTRAN implementation by Thomas Vojta, vojta@mst.edu
! built on a module found at www.fortran.com
!
!
! History:
!        v0.9     Dec 11, 2010    first implementation
!        V0.91    Dec 11, 2010    inlined internal function for the SR component
!        v0.92    Dec 13, 2010    extra shuffle of seed in kissinit
!        v0.93    Aug 13, 2012    changed integer representation test to avoid data statements
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      FUNCTION rkiss05()
      implicit none

      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers
      real(r8b),parameter    :: am=4.656612873077392578d-10       ! multiplier 1/2^31

      real(r8b)             :: rkiss05
      integer(i4b)          :: kiss
      integer(i4b)          :: x,y,z,w              ! working variables for the four generators
      common /kisscom/x,y,z,w

      x = 69069 * x + 1327217885
      y= ieor (y, ishft (y, 13)); y= ieor (y, ishft (y, -17)); y= ieor (y, ishft (y, 5))
      z = 18000 * iand (z, 65535) + ishft (z, - 16)
      w = 30903 * iand (w, 65535) + ishft (w, - 16)
      kiss = ishft(x + y + ishft (z, 16) + w , -1)
      rkiss05=kiss*am
      END FUNCTION rkiss05


      SUBROUTINE kissinit(iinit)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter     :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      integer(i4b) idum,ia,im,iq,ir,iinit
      integer(i4b) k,x,y,z,w,c1,c2,c3,c4
      real(r8b)    rkiss05,rdum
      parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
      common /kisscom/x,y,z,w

      !!! Test integer representation !!!
      c1=-8
      c1=ishftc(c1,-3)
!     print *,c1
      if (c1.ne.536870911) then
         print *,'Nonstandard integer representation. Stoped.'
         stop
      endif

      idum=iinit
      idum= abs(1099087573 * idum)               ! 32-bit LCG to shuffle seeds
      if (idum.eq.0) idum=1
      if (idum.ge.IM) idum=IM-1

      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         x=idum+1
      else
         x=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         y=idum+1
      else
         y=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         z=idum+1
      else
         z=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         w=idum+1
      else
         w=idum
      endif

      rdum=rkiss05()

      return
      end subroutine kissinit



 FUNCTION erfcc(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates complementary error function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      real(r8b)   :: erfcc,x
      real(r8b)   :: t,z
      z=abs(x)
      t=1.D0/(1.D0+0.5D0*z)
      erfcc=t*exp(-z*z-1.26551223D0+t*(1.00002368D0+t*(.37409196D0+t*&
     &(.09678418D0+t*(-.18628806D0+t*(.27886807D0+t*(-1.13520398D0+t*&
     &(1.48851587D0+t*(-.82215223D0+t*.17087277D0)))))))))
      if (x.lt.0.D0) erfcc=2.D0-erfcc
      return
      END  