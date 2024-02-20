# brain_FBM

The code in this project simulates reflected fractional Brownian motion (FBM) in complex two-dimensional and three-dimensional geometries. It was developed to simulate the spatial distribution of serotonergic neuros in vertebrate brains, see  [Front. Comput. Neurosci. 14, 56 (2020)](https://doi.org/10.3389/fncom.2020.00056) and [Front. Comput. Neurosci. 17, 1189853 (2023)](https://doi.org/10.3389/fncom.2023.1189853). It should be easy to adapt to other problems by modifying the files that describe the geometry (see below).

Each serotonergic axon is desribed as distinct discrete-time FBM trajectory, characterized by the decay exponent gamma. (gamma is related to the Hurst exponent via gamma = 2-2H). The increments (steps) of the FBM process, i.e., the fractional Gaussian noise, are created using the effective Fourier-filtering technique, allowing one to simulate long trajectories up to 2^25 (about 34 million) time steps. The code uses Ooura's FFT package, https://github.com/biotrump/OouraFFT, to perform the Fourier transformations. Conditional compilation using preprocessor directives is used to allow the same code to run serially (non-parallel) or in parallel using MPI. For code simulating the same reflected FBM process in simple geometries (interval, square, disk, sphere, etc) see [https://github.com/thomasvojta/reflected_FBM](https://github.com/thomasvojta/reflected_FBM) 

To achieve good performance, a fast way of checking whether a given point is inside or outside the given domain is crucial. In two dimensions, we therefore store the geometry data in a special format. The positions of the boundary are discretized, leading to integer coordinates. The outer boundary of the two-dimensioanl domain is stored in the file "outline_lr.dat". It contains an ordered list that specifies, for each y-coordinate, the minimum and maximum allowed x-values (i.e., the left and right boundaries). The format can handle convex domains only. Non-convex regions are treated as holes which are stored separately in the same format in the files "holeXX_lr.dat" where XX is the number of the hole. Such holes can also represent forbidden regions in the interior of the domain. The same principle is applied in 3 dimensions. The domain is first divided into two-dimensional crosssections. Each crosssection is stored as a 2d domain as described above. The outlines and holes are stored in the files shapXX_YY.dat and holeXX_YY.dat. Here YY counts the crosssections and XX numbers the outlines and holes in each crossection. The three-dimensional simulation than interpolates between the neighboring crosssections of a given point.

The main output of the simulations is the spatial probability density of the FBM walkers which models the spatial density of the serotonergic axons in the brain.

If you use this code, or any derivative work in academic work, please cite the following publications:
- Skirmantas Janušonis, Nils Detering, Ralf Metzler, and Thomas Vojta, Serotonergic Axons as Fractional Brownian Motion Paths: Insights Into the Self-Organization of Regional Densities, Front. Comput. Neurosci. 14, 56 (2020), [https://doi.org/10.3389/fncom.2020.00056](https://doi.org/10.3389/fncom.2020.00056)
- Skirmantas Janušonis, Justin H. Haiman, Ralf Metzler, and Thomas Vojta, Predicting the distribution of serotonergic axons: a supercomputing simulation of reflected fractional Brownian motion in a 3D-mouse brain model, Front. Comput. Neurosci. 17, 1189853 (2023), [https://doi.org/10.3389/fncom.2023.1189853](https://doi.org/10.3389/fncom.2023.1189853)

