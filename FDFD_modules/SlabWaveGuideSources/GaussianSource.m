
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
[TH,R] = cart2pol(X,Y);
[X,Y] = pol2cart(TH+theta,R);
nxs = round(0.3*Nx); %horizontal position of beam center
X = X - X(nxs,1);
fsrc = exp(-(2*X/w).^2).*exp(-1i*k0*n*Y);