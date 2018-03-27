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
sign = 1;
M = prod(N);  % total number of cells in domain
%take advantage of linear indexing

ind_cur = 1:M;  % indices of current points
ind_cur = ind_cur(:);

ind_adj_x = 1:M;  % indices of adjacent (previous or next) points in the w-direction
ind_adj_x = reshape(ind_adj_x, N);

%% shift the row indices
Dws = -4*sign*speye(M); %M fully determines the matricial size;
%% usual five point stencil
for w = ['x', 'y']
    ind_adj_r0x = circshift(ind_adj_x, -sign * ('xy' == w));
    ind_adj_l0x = circshift(ind_adj_x, sign * ('xy' == w));

    ind_adj_r0x = ind_adj_r0x(:);
    ind_adj_l0x = ind_adj_l0x(:);

    %% conver the offdiagonals into a linear index
    linear_ind_r0x = sub2ind([M M], ind_cur, ind_adj_r0x); 
    linear_ind_l0x = sub2ind([M M], ind_cur, ind_adj_l0x); 

    Dws(linear_ind_r0x) = sign;
    Dws(linear_ind_l0x) = sign;


end

%% rotated stencil (45 degrees)
Dws_r = -4*sign*speye(M,M); %M fully determines the matricial size;

ind_adj_rtx = circshift(ind_adj_x, sign * [1, 1]);
ind_adj_rbx = circshift(ind_adj_x, sign * [1, -1]);
ind_adj_ltx = circshift(ind_adj_x, sign * [-1, 1]);
ind_adj_lbx = circshift(ind_adj_x, sign * [-1, -1]);

ind_adj_rtx = ind_adj_rtx(:);
ind_adj_ltx = ind_adj_ltx(:);
ind_adj_rbx = ind_adj_rbx(:);
ind_adj_lbx = ind_adj_lbx(:);

%% conver the offdiagonals into a linear index
linear_ind_rtx = sub2ind([M M], ind_cur, ind_adj_rtx); 
linear_ind_ltx = sub2ind([M M], ind_cur, ind_adj_ltx); 
linear_ind_rbx = sub2ind([M M], ind_cur, ind_adj_rbx); 
linear_ind_lbx = sub2ind([M M], ind_cur, ind_adj_lbx); 

Dws_r(linear_ind_rtx) = sign;
Dws_r(linear_ind_ltx) = sign;
Dws_r(linear_ind_rbx) = sign;
Dws_r(linear_ind_lbx) = sign;


Dws = (1/dw^2)*Dws;
Dws_r = (1/(2*dw^2))*Dws_r;
figure;
spy(Dws_r)
figure;
spy(Dws)
a = 0.5;
Laplacian = a*Dws+(1-a)*Dws_r;

constant = 1000
A = Laplacian + constant*speye(N.^2);
% constant has to be tuned

b = zeros(length(Dws),1);
b(4750) = 1i*2

x1 = A\b;
E = reshape(x1, N(1), N(2));
visabs(E, [-1,1], [-1,1]);


