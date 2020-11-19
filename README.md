# Anisotropic-dispersive-2D-Kuramoto-Sivashinsky-Equation

An .exe, datafile and matlab script to perform and analyse numerical simulations.

Brief explanation of variables in the datafile (including one you should not change):

M, N are the truncations of the Fourier series, giving the size of the FFT matrix; these work best as powers of 2

dt is the time step, see "Linearly implicit schemes for multi-dimensional Kuramoto–Sivashinsky type equations arising in falling film flows" by Akrivis et al. for a convergence study (we use the second order BDF method). Larger dispersion means smaller dt is required

Row 3 should not be changed, alpha does not correspond to the alpha in the paper!

The code loops over the streamwise length L, and dispersion parameter delta, in rows 4 and 5 you can determine the starting values and number of increments. The matlab script is able to plot the results from the output data files one after another.

a corresponds to alpha in the paper, it controls the aspect ratio of the 2D domain, the transverse length is L^a, and with a=1 the simulations are on square domains.

NB: These simulations take a long time and output large data files!

NB2: Even though the sidebar says these codes are 100% Matlab, they are not, only the data handling and plotting tool "plots_mgd2DKSE.m" should be run in Matlab. Before this you must run the .exe file (how you run this is dependent on your computer), to solve the 2D KSE and generate data.
