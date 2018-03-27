
%% Gaussian Beam Source
Nx = 100;
Ny = 100;
xa = 1:Nx;
ya = 1:Nx;
theta = pi/2;
w = 2;
k0 = 2;
n = 1;

[Y,X] = meshgrid(ya,xa);
fsrc = exp(-1i*k0*n*Y);
omega = 1; 
c = 1;
k = [1,0];
J = PlaneWave(N, k, omega, c)
figure;
imagesc(abs(J));
