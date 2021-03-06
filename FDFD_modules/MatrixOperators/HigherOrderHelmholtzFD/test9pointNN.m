%% derivative operators with dirichlet boundary conditions
close all;
clear; 

Nx = 100; 
Ny = 100; 
N = [Nx, Ny]; dx = 0.01; dy = dx; dz = dx;
dL = [dx dy];
%% Solver Code Begins Here
w = 'x'; s = 'b';
%w = 'x', 'y', or 'z'
%s = 'b' or 'w'
% dL: [dx dy dz] for 3D; [dx dy] for 2D
% N: [Nx Ny Nz] for 3D; [Nx Ny] for 2D


dw = dL('xyz' == w);  % one of dx, dy, dz, all are numbers...no info about which one was selected
sign = 1;  % +1 for s=='f'; -1 for s=='b'
M = prod(N);  % total number of cells in domain
%take advantage of linear indexing

ind_cur = 1:M;  % indices of current points
ind_cur = ind_cur(:);

ind_adj_x = 1:M;  % indices of adjacent (previous or next) points in the w-direction
ind_adj_x = reshape(ind_adj_x, N);

%% shift the row indices
Dws = -(10/3)*sign*speye(M); %M fully determines the matricial size;
for w = ['x', 'y']
    ind_adj_r0x = circshift(ind_adj_x, -sign * ('xy' == w));
    ind_adj_l0x = circshift(ind_adj_x, sign * ('xy' == w));
    ind_adj_rtx = circshift(ind_adj_x, sign * [1, 1]);
    ind_adj_rbx = circshift(ind_adj_x, sign * [1, -1]);
    ind_adj_ltx = circshift(ind_adj_x, sign * [-1, 1]);
    ind_adj_lbx = circshift(ind_adj_x, sign * [-1, -1]);

    ind_adj_r0x = ind_adj_r0x(:);
    ind_adj_l0x = ind_adj_l0x(:);
    ind_adj_rtx = ind_adj_rtx(:);
    ind_adj_ltx = ind_adj_ltx(:);
    ind_adj_rbx = ind_adj_rbx(:);
    ind_adj_lbx = ind_adj_lbx(:);

    %% conver the offdiagonals into a linear index
    linear_ind_r0x = sub2ind([M M], ind_cur, ind_adj_r0x); 
    linear_ind_l0x = sub2ind([M M], ind_cur, ind_adj_l0x); 
    linear_ind_rtx = sub2ind([M M], ind_cur, ind_adj_rtx); 
    linear_ind_ltx = sub2ind([M M], ind_cur, ind_adj_ltx); 
    linear_ind_rbx = sub2ind([M M], ind_cur, ind_adj_rbx); 
    linear_ind_lbx = sub2ind([M M], ind_cur, ind_adj_lbx); 

    Dws(linear_ind_r0x) = (2/3)*sign;
    Dws(linear_ind_l0x) = (2/3)*sign;
    Dws(linear_ind_rtx) = (1/6)*sign;
    Dws(linear_ind_ltx) = (1/6)*sign;
    Dws(linear_ind_rbx) = (1/6)*sign;
    Dws(linear_ind_lbx) = (1/6)*sign;

end
Dws = (1/dw^2)*Dws;
%spy(Dws)
constant = 90;
omega = constant;

A = Dws + omega^2*speye(N.^2);
% constant has to be tuned

b = zeros(length(Dws),1);
b(4750) = 1i*omega*1;

x1 = A\b;
figure; 
E = reshape(x1, N(1), N(2));
visabs(E, [-1,1], [-1,1]);

%% compare with the usual 5 point formulation
Dxf = createDws_dense('x', 'f', dL, N); 
Dyf = createDws_dense('y', 'f', dL, N);
Dyb = createDws_dense('y', 'b', dL, N); 
Dxb = createDws_dense('x', 'b', dL, N); 
A2 = Dxf*Dxb+Dyf*Dyb + omega^2*speye(N.^2);

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

