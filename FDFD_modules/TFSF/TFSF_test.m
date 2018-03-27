close all
clear

%% specify a simple FDFD
c0 = 3*10^8;
L0 = 1e-6;  % length unit: microns
wvlen = 5.0;  % wavelength in L0
xrange = [-5 5];  % x boundaries in L0
yrange = [-5 5];  % y boundaries in L0
Npml = 0*[15 15];  % [Nx_pml Ny_pml]
N = [200,200];

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
%% TFSF Partition
Q = ones(N); %0 indexes the total field
Q(25:75, 25:75) = 0; %1 indexes the scattered field, so we isolate only scattered field
Q(45:55, 45:55) = 1;
% convert into diagonal matrix sparse(
M = prod(N);
Qr = diag(sparse(Q(:)));

%% Set up the permittivity.
eps_r = ones(N);

%% Set up the 1agnetic current source density.
%must be entirely inside the TFSF interface
ind_src = [50,50];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
xa = 1:N(1);
ya = 1:N(2);
[Y,X] = meshgrid(ya,xa);
h = dL(1);
k0 = 2*pi/wvlen;
Jz = exp(-1i*(h *k0*(X)));

%%solve for the source field in vacuum
[Ez, Hx, Hy, A0, omega0,b] = solveTM(wvlen, xrange, yrange, eps_r, Jz, Npml);
fsrc = A0\b; %vacuum field distribution of the point source


Fsrc = reshape(fsrc, Nx,Ny);
figure; 
visreal(Fsrc, xrange, yrange);

%% ADD SCATTERER
eps_r(20:60, 20:60) = 1;
[Eza, Hxa, Hya, A, omega,b] = solveTM(wvlen, xrange, yrange, eps_r, Jz, Npml);

bprime = (Qr*A - A*Qr)*fsrc;

%% solve TFSF problem
x = A\bprime; %% the tfsf problem isolates ONLY  the scattered fields...?

e = reshape(x, Nx,Ny);
visreal(1i*e, xrange, yrange);

figure;
visreal(1i*Eza, xrange, yrange);

E_TF = TFSF(Qr,A,Jz,L0,wvlen,xrange, yrange, Npml);
figure; 
visreal(1i*E_TF, xrange, yrange);

