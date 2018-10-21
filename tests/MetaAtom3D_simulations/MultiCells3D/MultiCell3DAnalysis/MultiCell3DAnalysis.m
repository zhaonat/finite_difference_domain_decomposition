
%% Parameter Set up
%% Parameter Set-up
L0 = 1e-6;  % length unit: microns
wvlen = 15.0;  % wavelength in L0
xrange = [-5 5];  % x boundaries in L0
yrange = [-5 5];  % y boundaries in L0
zrange = [-5 5];
Nx = 20; Ny = 20; Nz = 20;
N = [Nx Ny Nz];  % [Nx Ny]
M = N(1)*N(2)*N(3);
Npml = [0 0 0];  % [Nx_pml Ny_pml]

%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  
% domain is expanded to include PML

%% Set up the permittivity.
eps_r = ones(N);

%% Set up the current source density.
Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = [5 5 5];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
JCurrentVector = [Mx; My; Mz];
%% Wonsoek's scalar parameter 1, -1, or 0
s = -1;
%=========================================

[Ex, Ey, Ez, A,b, omega, c0, Tepsuper, Tmusuper] = ...
    solve3D_EigenEngine(L0, wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml,s);

%MultiCellReorder(A, b,xCells, yCells, N)
%% parameters for a multicellreordering
xCells = 2; yCells = 1; zCells = 3;
Grid = zeros(Nx, Ny, Nz);
counter = 0;
%% Solver Code Starts Here

    [permutedIndices, boundaryCells, interiorCells, hpart, vpart] = ...
    IndexPermutation3DTwoDCubicLattice(N,xCells, yCells);

   
    %% Execute a Symmetry Preserving Column Row Permutation Combination
    Q = speye(3*N(1)*N(2)*N(3));
    %[linesF, columnsF, valuesF] = find(Q);
    
    for i = 1:3*N(1)*N(2)*N(3)
       Q(i,i) = 0;
       indexshift = permutedIndices(i);
%        linesF = [linesF; i]; THIS PRODUCES VERY INTERESTING MATRICES
%        columnsF = [columnsF; indexshift];
%        valuesF = [valuesF; 1];
       Q(i,indexshift) = 1;
    end
    %Q = sparse(linesF, columnsF, valuesF);

    %%test permutation matrix should permute index order to permuted indices

    %% Transform Equations
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    SymB = Q*b;
    
   %% Execute IterativeFourBlockSchur
   
   [Aschur, bmod] = ...
    FourblockSchur(SymA,SymB, hpart,vpart);

figure;
spy(SymA);

figure;
plot(abs(Aschur\bmod));
 
