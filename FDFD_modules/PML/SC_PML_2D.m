%% generate scpml

function [sxf_grid, syf_grid, sxb_grid, syb_grid] = ...
    SC_PML_2D(xrange, yrange, omega, eps_0, mu_0,...
    N, Npml)
    M = prod(N);
    Nx_pml = Npml(1);
    Ny_pml = Npml(2);
    Nx = N(1);
    Ny = N(2);
    sxf = create_sfactor_mine(xrange,'f',omega,eps_0,mu_0,Nx,Nx_pml);
    syf = create_sfactor_mine(yrange,'f', omega,eps_0,mu_0,Ny,Ny_pml);
    sxb = create_sfactor_mine(xrange, 'b', omega,eps_0,mu_0, Nx, Nx_pml);
    syb = create_sfactor_mine(yrange,'b', omega,eps_0,mu_0,Ny,Ny_pml);

    % now we create the matrix (i.e. repeat sxf Ny times repeat Syf Nx times)
    [Sxf, Syf] = ndgrid(sxf, syf);
    [Sxb, Syb] = ndgrid(sxb, syb);
    sxf_grid = Sxf; syf_grid = Syf;
    sxb_grid = Sxb; syb_grid = Syb;
    %Sxf(:) converts from n x n t0 n^2 x 1
    Sxf=spdiags(Sxf(:),0,M,M);
    Sxb=spdiags(Sxb(:),0,M,M);
    Syf=spdiags(Syf(:),0,M,M);
    Syb=spdiags(Syb(:),0,M,M);

end