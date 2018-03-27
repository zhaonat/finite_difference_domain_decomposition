%% TEST higher order FD
close all;
clear; 


Nx = 100; 
Ny = 100; 
N = [Nx, Ny]; dx = 0.01; dy = dx; dz = dx;
dL = [dx dy];
omega = 100;

%% compare with the usual 5 point formulation
Dxf = createDws_dense('x', 'f', dL, N); 
Dyf = createDws_dense('y', 'f', dL, N);
Dyb = createDws_dense('y', 'b', dL, N); 
Dxb = createDws_dense('x', 'b', dL, N); 
A2 = Dxf*Dxb+Dyf*Dyb + omega^2*speye(N.^2);
b = zeros(length(A2),1);
b(4750) = 1i*omega*1;
x2 = A2\b;
figure; 
E2 = reshape(x2, N(1), N(2));
visabs(E2, [-1,1], [-1,1]);
k = omega; h = dL(1);
A_h = Helmholtz4thOrder(omega, N,dL);
b4 = (2/3-(k*h)^2)/12*b;
x3 = A_h\b4;

figure; 
E3 = reshape(x3, N(1), N(2));
visabs(E3, [-1,1], [-1,1]);

A_h6 = Helmholtz6thOrder(omega, N,dL);
x3 = A_h6\b;

figure; 
E3 = reshape(x3, N(1), N(2));
visabs(E3, [-1,1], [-1,1]);

