      PROGRAM fbm_xy_lr7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   time-discrete reflected fractional Brownian motion
!   on 2D shape read from file
!
!   27 Apr 2019      fbm_xy_lr1       first version, creates individual trajectories
!   27 Apr 2019      fbm_xy_lr2       measure 2d density averaged over several trajectories
!   27 Apr 2019      fbm_xy_lr3       define function for inside/outside
!   28 Apr 2019      fbm_xy_lr4       complete implementation of parallization 
!   30 Apr 2019      fbm_xy_lr5       add uncorrelated disorder
!   30 Apr 2019      fbm_xy_lr6       three obstacles
!   15 Jun 2019      fbm_xy_lr7       variable number of holes, reads ymax, ymin   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preprocessor directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#define PARALLEL
#define VERSION 'fbm_xy_lr7'

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
      integer(i4b), parameter     :: NCONF=960                    ! number of walkers
      real(r8b), parameter        :: GAMMA = 0.4D0              ! FBM correlation exponent 
      logical, parameter          :: FRACTIONAL=.true.           ! true: FBM, false: BM

      integer(i4b),parameter      :: NHOLE=6                    ! number of holes in the shape 
      
      real(r8b),parameter         :: X0=110.D0, Y0=100.D0          ! starting point
      logical, parameter          :: RANDOMSTART=.true.          ! if true start from random position

      real(r8b), parameter        :: STEPSIG=0.1D0                ! sigma of individual step 
      character(3), parameter     :: STEPDIS='GAU'                 ! Gaussian = GAU, binary = BIN, box = BOX                   

      real(r8b), parameter        :: GRIDXMIN=4.D0, GRIDXMAX=386.D0   ! grid dimensions for measuring density
      real(r8b), parameter        :: GRIDYMIN=14.D0, GRIDYMAX=342.D0
      real(r8b), parameter        :: CELLSIZE=1.D0                      ! size of grid cell  
       
      
      logical, parameter          :: WRITETRAJEC = .false.        ! write individual trajectories
      logical, parameter          :: WRITEDISTRIB = .true.       ! write density distribution
      
      integer(i4b), parameter     :: IRINIT=1                    ! random number seed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      integer(i4b), parameter     :: NGX=NINT((GRIDXMAX-GRIDXMIN)/CELLSIZE)          ! number of grid cells in x direction
      integer(i4b), parameter     :: NGY=NINT((GRIDYMAX-GRIDYMIN)/CELLSIZE)          ! number of grid cells in y direction
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(r8b)                 :: xx(0:NT), yy(0:NT)               ! walker coordinates
      real(r8b)                 :: xix(1:NT), xiy(1:NT)             ! increments
      real(r8b)                 :: rad2                             ! distance of walker from origin
      integer(i4b)              :: nintyy                           ! nint(yy(it)) 
      
      integer(i4b)              ::ymin(0:NHOLE), ymax(0:NHOLE) 
      integer(i4b), allocatable :: xleft(:,:), xright(:,:)
      
      logical                :: outside
      logical                :: inhole(NHOLE)
      
      real(r8b)              :: confxx(1:NT)                     ! average of xx after each time step  
      real(r8b)              :: conf2xx(1:NT)                    
      real(r8b)              :: sumxx(1:NT),sum2xx(1:NT)         ! sums over machines in MPI version
      real(r8b)              :: auxxx(1:NT),aux2xx(1:NT) 

      integer(i4b)           :: iconf, it, iy, ihole             ! configuration, time, and hole counters   
      integer(i4b)           :: totconf                         ! actual number of confs
      
      real(r8b)              :: confdens(NGX,NGY)                ! density array
      real(r8b)              :: normdens, lognormdens             ! normalized density 
      real(r8b)              :: cellx,celly                      ! center coordinates of each cell
      integer(i4b)           :: igx,igy                          ! grid cell counter
      real(r8b)              :: sumdens(NGX,NGY)                ! sums over machines in MPI version
      real(r8b)              :: auxdens(NGX,NGY)                ! 
      
      external               :: kissinit 
      real(r8b),external     :: rkiss05,gkiss05, erfcc 

      character(14)          :: outlinefile = "outline_lr.dat"
      character(13)          :: holefile = "holeXX_lr.dat"
     
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

! Read shape information 
      open(8,file=outlinefile,status='old')
      rewind(8)
      read(8,*) cdummy, ymin(0)
      read(8,*) cdummy, ymax(0)
!      print *, ymin(0),ymax(0)
      allocate( xleft(0:NHOLE,ymin(0):ymax(0)),xright(0:NHOLE,ymin(0):ymax(0)) )
      do iy=ymin(0),ymax(0)
         read (8,*) dummy, xleft(0,iy),xright(0,iy)
      enddo
      close(8) 
      
      do ihole=1,NHOLE
        write(holefile(5:6),'(I2.2)') ihole 
        open(8,file=holefile,status='old')
        rewind(8)
        read(8,*) cdummy, ymin(ihole)
        read(8,*) cdummy, ymax(ihole)
!        print *, ymin(ihole),ymax(ihole)
          do iy=ymin(ihole),ymax(ihole)
            read (8,*) dummy, xleft(ihole,iy),xright(ihole,iy)
          enddo
        close(8)
      enddo
      
      confxx(:)=0.D0 
      conf2xx(:)=0.D0 

      confdens(:,:)=0.D0
      
! Loop over disorder configurations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      disorder_loop: do iconf=myid+1,totconf,numprocs
!      if (myid==0) print *, 'dis. conf.', iconf
#else
      disorder_loop: do iconf=1,totconf
!         print *, 'dis. conf.', iconf
#endif 

        call gkissinit(IRINIT+iconf-1)

        if (FRACTIONAL) then
          call corvec(xix,NT,M,GAMMA)                         ! create x-increments (correlated Gaussian random numbers)
          call corvec(xiy,NT,M,GAMMA)                         ! create y-increments (correlated Gaussian random numbers)
        else
          do it=1,NT  
            xix(it)=gkiss05()                ! create uncorrelated random step
            xiy(it)=gkiss05()                ! create uncorrelated random step
          enddo 
        endif  
        
                
        if (STEPDIS.eq.'BOX') then
          xix(:) = 1 - (0.5D0*erfcc(xix(:)/sqrt(2.0D0)))                       ! map onto unit inteval
          xix(:)=xix(:)-0.5D0                                                  ! center around 0
          xiy(:) = 1 - (0.5D0*erfcc(xiy(:)/sqrt(2.0D0)))                       ! map onto unit inteval
          xiy(:)=xiy(:)-0.5D0                                                  ! center around 0
        endif
        
        if (STEPDIS.eq.'BIN') then 
          xix(:)=sign(1.D0,xix(:))
          xiy(:)=sign(1.D0,xiy(:))
        endif  
        
        xix(:)=xix(:)*STEPSIG
        xiy(:)=xiy(:)*STEPSIG
        
! Time loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           call initialize_point 
           do it=1, NT/2
              xx(it)=xx(it-1)+xix(it)
              yy(it)=yy(it-1)+xiy(it)
      
              if (pointnotvalid(xx(it),yy(it))) then
                 xx(it)=xx(it-1)
                 yy(it)=yy(it-1)
              end if
              
              igx= 1+int(xx(it)-GRIDXMIN)/CELLSIZE
              igy= 1+int(yy(it)-GRIDYMIN)/CELLSIZE
              confdens(igx,igy)=confdens(igx,igy)+1.D0
           end do  

        if (WRITETRAJEC) then   
          write(trafile(4:9),'(I6.6)') IRINIT+iconf-1       
          open(2,file=trafile,status='replace')
          write(2,*) 'Program ', VERSION
          write(2,*) 'trajectory of reflected FBM on 2d shape'
          write(2,*) 'step distribution ',STEPDIS
          write(2,*) 'step size ',STEPSIG
          write(2,*) 'NT= ', NT, ' number of steps used:', NT/2
          write(2,*) 'FRACTIONAL=', FRACTIONAL
          write(2,*) 'GAMMMA=', GAMMA
          write(2,*) 'IRINIT=',IRINIT
          write(2,*) 'RNG seed=', IRINIT+iconf-1
          write (2,*)'=================================='
          write(2,*) '   time         x         y'
          do it=0, NT/2
            Write(2,'(1X,I7,6(2X,E13.6))')  it, xx(it), yy(it)
          enddo 
          close(2) 
        endif  

      end do disorder_loop      ! of do incof=1,NCONF
      
! Now collect and communicate data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
     if (myid.ne.0) then                                                  ! Send data
         call MPI_SEND(confdens,NGX*NGY,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
      else
         sumdens(:,:)=confdens(:,:)
         do id=1,numprocs-1                                                   ! Receive data
            call MPI_RECV(auxdens,NGX*NGY,MPI_DOUBLE_PRECISION,id,1,MPI_COMM_WORLD,status,ierr)
            sumdens(:,:)=sumdens(:,:)+auxdens(:,:) 
         enddo
      endif        
#else          
      sumdens(:,:)=confdens(:,:)
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
          write(2,*) 'density distribution of reflected FBM on 2d shape'
          write(2,*) 'step distribution ',STEPDIS
          write(2,*) 'step size ',STEPSIG
          write(2,*) 'NT= ', NT,' number of steps used:', NT/2  
          write(2,*) 'FRACTIONAL=', FRACTIONAL
          write(2,*) 'GAMMMA=', GAMMA
          write(2,*) 'Number of walkers', totconf
          write(2,*) 'IRINIT=',IRINIT
          write (2,*)'=================================='
          write(2,*) '   igx  igy    x      y     distrib    ln(distrib)'
          do igx=1,NGX
          do igy=1,NGY
            cellx= GRIDXMIN + (1.D0*igx -0.5D0)*CELLSIZE
            celly= GRIDyMIN + (1.D0*igy -0.5D0)*CELLSIZE
            normdens=sumdens(igx,igy)/((1.D0*totconf)*(NT/2)*CELLSIZE**2)
            if (normdens.gt.0.D0) then
               lognormdens=log(normdens)
            else
               lognormdens=-100
            endif   
            Write(2,'(1X,I4,1x,I4,4(1X,E13.6))')  igx, igy, cellx, celly, normdens, lognormdens
          enddo 
          enddo
          close(2)             
      endif
#ifdef PARALLEL
      endif ! of if (myid==0)
#endif        

      deallocate(xright,xleft)
      
#ifdef PARALLEL
      call MPI_FINALIZE(ierr)
#endif 
      stop      
      
! Now the internal subroutines !!!!!!!!!!!!!!!!!!!!!
      contains
            
      function pointnotvalid(xn,yn)
      real(r8b) :: xn,yn
      logical   :: pointnotvalid
      nintyy=nint(yn)
      
      if ( nintyy < ymin(0) .or.  nintyy > ymax(0) ) then
         pointnotvalid=.true.
      else    
        outside =  xn < xleft(0,nintyy) .or. xn > xright(0,nintyy)
        do ihole=1,NHOLE
          inhole(ihole) = nintyy >=ymin(ihole) .and. nintyy<=ymax(ihole) .and. xn>xleft(ihole,nintyy) .and. xn<xright(ihole,nintyy)  
        enddo  
        pointnotvalid=outside 
        do ihole=1,NHOLE
           pointnotvalid=pointnotvalid .or. inhole(ihole) 
        enddo
      endif  
      end function pointnotvalid

      subroutine initialize_point
      if (RANDOMSTART) then            
         xx(0)=0.D0
         yy(0)=0.D0  
         do while (pointnotvalid(xx(0),yy(0)))
           xx(0)=GRIDXMIN+rkiss05()*(GRIDXMAX-GRIDXMIN) 
           yy(0)=GRIDYMIN+rkiss05()*(GRIDYMAX-GRIDYMIN) 
         enddo
      else
         xx(0)=X0
         yy(0)=Y0  
      endif   

      end subroutine initialize_point
            
      END PROGRAM fbm_xy_lr7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine CORVEC(xr,Ns,M)
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

      SUBROUTINE corvec(xr,Ns,M,gam)
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
            
      real(r8b), external    :: gkiss05,erfcc
      external               :: rdft                   ! from Ooura's FFT package
      real(r8b),external     :: fbmcorfunc 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (Ns.ne.2**M) STOP 'Size indices do not match'
! Calculate correlation function 
      do is=0,Ns-1 
         cr(is)= fbmcorfunc(is,Ns,gam) 
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
      FUNCTION fbmcorfunc(i,N,gam)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      real(r8b)              :: corfunc,fbmcorfunc
      real(r8b)              :: gam
      integer(i4b)           :: i,N
      integer(i4b)           :: dist

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

      fbmcorfunc=corfunc
      
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