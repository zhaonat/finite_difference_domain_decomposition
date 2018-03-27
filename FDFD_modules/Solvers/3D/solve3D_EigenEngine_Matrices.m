
function [A, b, omega, c0, TepsSuper, Tmu, ...
    Sxf, Syf, Szf, Sxb, Syb, Szb] = ...
    solve3D_EigenEngine_Matrices(L0, wvlen, xrange, yrange, zrange, eps_r,...
    JCurrentVector, Npml, s)

    %% Set up the domain parameters.
    %normal SI parameters
    eps_0 = 8.85*10^-12*L0;
    mu_0 = 4*pi*10^-7*L0; 
    c0 = 1/sqrt(eps_0*mu_0);  % speed of light in 
    N = size(eps_r);  
    omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

    %% Set up the permittivity and permeability in the domain.
    % bwdmean does nearest neighbor averaging (smoothes out stuff)
    eps_x = bwdmean_w(eps_0 * eps_r, 'y');  % average eps for eps_x
    eps_y = bwdmean_w(eps_0 * eps_r, 'x');  % average eps for eps_y
    eps_z = bwdmean_w(eps_0 * eps_r, 'z');

    %% Set up number of cells
    %the wavelength is much larger than the dimensions of the system...
    xmin = xrange(1); xmax = xrange(2);
    ymin = yrange(1); ymax = yrange(2);
    zmin = zrange(1); zmax = zrange(2);
    Nx = N(1); dx = abs(xmax - xmin)/Nx;
    Ny = N(2); dy = abs(ymax - ymin)/Ny;
    Nz = N(3); dz = abs(zmax - zmin)/Nz;
    % Nz = 1; dz = 1; 2D solving only
    M = prod([Nx, Ny, Nz]); %total number of cells
    Mx = JCurrentVector(1:Nx,:,:);
    My = JCurrentVector(Nx+1:2*Nx,:,:);
    Mz = JCurrentVector(2*Nx+1: 3*Nx,:,:);
    
    %% Create the dielectric and permeability arrays (ex, ey, muz)
    %create a diagonal block matrix of ep and mu...
    epxList = reshape(eps_x,M,1);
    epyList = reshape(eps_y,M,1);
    epzList = reshape(eps_z,M,1);
    Teps_x = spdiags(epxList,0,M,M); % creates an MxM matrix, which is the correct size,
    %the M entries in epsList is put on the diagonals
    Teps_y = spdiags(epyList,0,M,M);
    Teps_z = spdiags(epzList, 0,M,M);
    Teps = blkdiag(Teps_x, Teps_y, Teps_z);
    Tmz = mu_0*speye(M); %in most cases, permeability is that of free-space

    %% Create SuperVector of the Dielectrics and Permeability
    TepsSuper = blkdiag(Teps_x,Teps_y,Teps_z);
    Tmu = blkdiag(Tmz, Tmz, Tmz);
    %% Create Current Source vector J 
    % dimension = M*1
    Mz = reshape(Mz, M, 1);
    My = reshape(My, M, 1);
    Mx = reshape(Mx, M, 1);
    Mz = sparse(Mz);

    %% create the derivative oeprators w/ PML

    N = [Nx, Ny, Nz];
    dL = [dx dy dz]; % Remember, everything must be in SI units beforehand

    Dxf = createDws_dense('x', 'f', dL, N); 
    Dxb = createDws_dense('x', 'b', dL, N); 
    Dyf = createDws_dense('y', 'f', dL, N);
    Dyb = createDws_dense('y', 'b', dL, N); 
    Dzf = createDws_dense('z', 'f', dL, N); 
    Dzb = createDws_dense('z', 'b', dL, N); 
    
    %% ADD IN PML
    sxf = create_sfactor_mine(xrange,'f',omega,eps_0,mu_0,Nx,Npml(1));
    sxb = create_sfactor_mine(xrange, 'b', omega,eps_0,mu_0, Nx, Npml(1));
    syf = create_sfactor_mine(yrange,'f', omega,eps_0,mu_0,Ny,Npml(2));
    syb = create_sfactor_mine(yrange,'b', omega,eps_0,mu_0,Ny,Npml(2));
    szf = create_sfactor_mine(zrange, 'f', omega,eps_0,mu_0, Nz, Npml(3));
    szb = create_sfactor_mine(zrange,'b', omega,eps_0,mu_0,Nz, Npml(3));

    [Sxf, Syf, Szf] = ndgrid(sxf, syf, szf);
    [Sxb, Syb, Szb] = ndgrid(sxb, syb, szb);
    
    %Sxf(:) converts from n x n t0 n^2 x 1
    Sxf=spdiags(Sxf(:),0,M,M);
    Sxb=spdiags(Sxb(:),0,M,M);
    Syf=spdiags(Syf(:),0,M,M);
    Syb=spdiags(Syb(:),0,M,M);
    Szf=spdiags(Szf(:),0,M,M);
    Szb=spdiags(Szb(:),0,M,M);
    
    Dxf = Sxf^-1*Dxf; 
    Dyf = Syf^-1*Dyf;
    Dyb = Syb^-1*Dyb; 
    Dxb = Sxb^-1*Dxb;
    Dzf = Szf^-1*Dzf;
    Dzb = Szf^-1*Dzb;
    
    %% CONSTRUCT SPECIAL AVERAGING TERM
    eps_xy = (circshift(eps_x, [1 0 0]) + eps_x)/2;
    Teps_xy = spdiags(reshape(eps_xy,M,1), 0, M, M);
    
    %% Construct Ch and Ce operators
    Ce = [sparse(M,M), -Dzf, Dyf; Dzf, sparse(M,M), -Dxf; -Dyf, Dxf, sparse(M,M)];
    Ch = [sparse(M,M), -Dzb, Dyb; Dzb, sparse(M,M), -Dxb; -Dyb, Dxb, sparse(M,M)];

    %% constrct the eE term
    %gradient(divergence)
%     GradDiv = [Dxf*Teps_xy^-1*Dxb*Tepx, Dxf*Teps_xy^-1*Dyb*Tepy, Dxf*Teps_xy^-1*Dzb*Tepz; ...
%         Dyf*Tepy^-1*Dxb*Tepx, Dyf*Tepy^-1*Dyb*Tepy, Dyf*Tepy^-1*Dzb*Tepz; ...
%         Dzf*Tepz^-1*Dxb*Tepx, Dzf*Tepz^-1*Dyb*Tepy, Dzf*Tepz^-1*Dzb*Tepz];
    D_D_j = [Dxf*Teps_xy^-1*Dxb,Dxf*Teps_xy^-1*Dyb,Dxf*Teps_xy^-1*Dzb;...
        Dyf*Teps_xy^-1*Dxb,Dyf*Teps_xy^-1*Dyb,Dyf*Teps_xy^-1*Dzb;...
        Dzf*Teps_xy^-1*Dxb,Dzf*Teps_xy^-1*Dyb,Dzf*Teps_xy^-1*Dzb];

    D_D_e = [Dxf*Teps_xy^-1*Dxb*Teps_x,Dxf*Teps_xy^-1*Dyb*Teps_y,Dxf*Teps_xy^-1*Dzb*Teps_z;...
        Dyf*Teps_xy^-1*Dxb*Teps_x,Dyf*Teps_xy^-1*Dyb*Teps_y,Dyf*Teps_xy^-1*Dzb*Teps_z;...
        Dzf*Teps_xy^-1*Dxb*Teps_x,Dzf*Teps_xy^-1*Dyb*Teps_y,Dzf*Teps_xy^-1*Dzb*Teps_z];

    J = [Mx; My; Mz];
    b = -1i * omega * Tmu * J  + s * 1i / omega * D_D_j * J;
    A = Ch * Ce - omega^2 * Tmu * Teps + s * D_D_e;

    %% construct the matrix b, everything is in 2D matrices
%     b = -1i*omega*J; bo = b;
%     JCorrection = (1i/omega) * (s*WAccelScal)*TepsSuper^-1*GradDiv*J;
%     b = b + JCorrection;
   
end