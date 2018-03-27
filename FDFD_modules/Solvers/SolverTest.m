%% Solver Testing
clear;
%% Set up the domain parameters.
c0 = 3*10^8;
L0 = 1e-6;  % length unit: microns
wvlen = 4.0*L0;  % wavelength in L0
xrange = [-5 5]*L0;  % x boundaries in L0
yrange = [-5 5]*L0;  % y boundaries in L0
N = [100 100];  % [Nx Ny]
Npml = [10 10];  % [Nx_pml Ny_pml]
mu0 = 4*pi*10^-7; mu_0 = mu0; mu = mu0;
eps0 = 8.85*10^-12; eps_0 = eps0;
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
eps_r = ones(N);

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = [25 25];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;


%% Run the solveTE Function
[Hz, Ex, Ey, A, omega,b, Sxf, Dxf, Dyf, sxf, syf] = solveTE(wvlen, xrange, yrange, eps_r, Mz, Npml);
visabs(Hz, xrange, yrange)
figure
visabs(Ex, xrange, yrange)
figure
visabs(Ey, xrange, yrange)
figure;
moviereal(Hz, xrange, yrange)