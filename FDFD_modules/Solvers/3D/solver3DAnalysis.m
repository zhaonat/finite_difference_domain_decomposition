%% Solver 3D Analysis
%% Parameter Set-up
L0 = 1e-6;  % length unit: microns
wvlen = 1.0;  % wavelength in L0
xrange = [-1 1];  % x boundaries in L0
yrange = [-1 1];  % y boundaries in L0
zrange = [-1 1];
Nx = 30; Ny = 30; Nz = 30;
N = [Nx Ny Nz];  % [Nx Ny]
M = N(1)*N(2)*N(3);
Npml = [0 0 0];  % [Nx_pml Ny_pml]

%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
eps_r = ones(N);

%% Set up the current source density.
Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = [Nx/2 Ny/2 Nz/2];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd

My(ind_src(1), ind_src(2), ind_src(3)) = 1;
JCurrentVector = [Mx; My; Mz];

%% SOLVER3D RAW CODE STARTS HERE

    %normal SI parameters
    eps_0 = 8.85*10^-12*L0;
    mu_0 = 4*pi*10^-7*L0; 
    eps0 = eps_0;  % vacuum permittivity
    mu0 = mu_0;  % vacuum permeability in
    c0 = 1/sqrt(eps_0*mu_0);  % speed of light in 
    N = size(eps_r);  % [Nx Ny] THIS IS THE POINT WHERE THE GRID SIZE IS DETERMINED
    omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

    %% Set up the permittivity and permeability in the domain.
    % bwdmean does nearest neighbor averaging (smoothes out stuff)
    eps_x = bwdmean_w(eps0 * eps_r, 'y');  % average eps for eps_x
    eps_y = bwdmean_w(eps0 * eps_r, 'x');  % average eps for eps_y
    eps_z = bwdmean_w(eps0 * eps_r, 'z');
    %these are fully dense matrices...

    %currently, eps_x and eps_y are ultra-dense, which isn't right...

    %% Set up number of cells
    %the wavelength is much larger than the dimensions of the system...
    xmin = xrange(1); xmax = xrange(2);
    ymin = yrange(1); ymax = yrange(2);
    zmin = zrange(1); zmax = zrange(2);
    Nx = N(1); dx = (xmax - xmin)/Nx;
    Ny = N(2); dy = (ymax - ymin)/Ny;
    Nz = N(3); dz = (zmax - zmin)/Nz;
    % Nz = 1; dz = 1; 2D solving only
    
    M = prod([Nx, Ny, Nz]); %total number of cells
    Mx = JCurrentVector(1:Nx,:,:); My = JCurrentVector(Nx+1:2*Nx,:,:); 
    Mz = JCurrentVector(2*Nx+1: 3*Nx,:,:);

    %% Create the dielectric and permeability arrays (ex, ey, muz)
    %create a diagonal block matrix of ep and mu...
    epxList = reshape(eps_x,M,1);
    epyList = reshape(eps_y,M,1);
    epzList = reshape(eps_z,M,1);
    Tepx = spdiags(epxList,0,M,M); % creates an MxM matrix, which is the correct size,
    %the M entries in epsList is put on the diagonals
    Tepy = spdiags(epyList,0,M,M);
    Tepz = spdiags(epzList, 0,M,M);
    Tmz = mu0*speye(M); %in most cases, permeability is that of free-space

    %% Create SuperVector of the Dielectrics and Permeability
    TepsSuper = blkdiag(Tepx,Tepy,Tepz);
    TmuSuper = blkdiag(Tmz, Tmz, Tmz);
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
    Dyf = createDws_dense('y', 'f', dL, N);
    Dyb = createDws_dense('y', 'b', dL, N); 
    Dxb = createDws_dense('x', 'b', dL, N); 
    Dzf = createDws_dense('z', 'f', dL, N); 
    Dzb = createDws_dense('z', 'b', dL, N); 

    %% Construct Ch and Ce operators
    Ce = [zeros(M) -Dzf Dyf; Dzf zeros(M) -Dxf; -Dyf Dxf zeros(M)];
    Ch = [zeros(M) -Dzb Dyb; Dzb zeros(M) -Dxb; -Dyb Dxb zeros(M)];

    %% Construct the matrix A
    A = Ch*TmuSuper^-1*Ce - omega^2*TepsSuper;
    
%     figure; 
%     spy(A); pause; 

    %% construct the matrix b, everything is in 2D
    J = [Mx; My; Mz];
    b = -1i*omega*J;
    solution = qmr(A,b,1e-8, length(b));
    solLength = length(solution);
    Ex = solution(1:solLength/3);
    Ey = solution(solLength/3+1:solLength*(2/3));
    Ez = solution(solLength*(2/3)+1: solLength);
    ExC = reshape(full(Ex), N(1), N(2), N(3));
    EyC = reshape(full(Ey), N(1), N(2), N(3));
    EzC = reshape(full(Ez), N(1), N(2), N(3));

%% Visualize Solution
figure;
EzC = (abs(EzC));
xslices = [10 20];
slice(log(EzC), xslices, xslices, xslices)
colormap b2r
colorbar
xlabel('x direction')
ylabel('y direction')
zlabel('z direction')
