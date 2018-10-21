close all
clear;

%% Parameter Set up
%% Parameter Set-up
L0 = 1e-6;  % length unit: microns
wvlen = 2.0;  % wavelength in L0

%% Create Data File Directory for this simulation
cd(strcat('D:\Nathan\Documents\StanfordYearOne\Fan Group\FDFDSchurComplement\MetaAtom3D\MultiCells3D\MultiCellIterations'));

dataDir = strcat('MultiCell3D 1D Chain with SingleSellSize = ',int2str(15),...
    'wvlen=',num2str(wvlen));
mkdir(dataDir);

cd(strcat('D:\Nathan\Documents\StanfordYearOne\Fan Group\FDFDSchurComplement\MetaAtom3D\MultiCells3D\MultiCellIterations\',dataDir));

iterations = [];
for i = 3:3
    Nx = 10*i; Ny = 10; Nz = 10;
    N = [Nx Ny Nz];  % [Nx Ny]
    M = N(1)*N(2)*N(3);
    Npml = [0 0 0];  % [Nx_pml Ny_pml]
    xrange = [-1 1];  % x boundaries in L0
    yrange = [-1 1];  % y boundaries in L0
    zrange = [-1 1];
    %[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  
    % domain is expanded to include PML

    %% Set up the permittivity.
    eps_r = ones(N);
    eps_r(2:8, 2:8, 2:8) = 12; eps_r(12:18, 2:8, 2:8) = 12;

    %% Set up the current source density.
    Mz = zeros(N); My = Mz; Mx = Mz;
    ind_src = [5 5 5];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
    JCurrentVector = [Mx; My; Mz];

    %% Wonsoek's scalar parameter 1, -1, or 0
    s = -1;
    %=========================================
    tic
    [A,b,Ao, bo, omega, c0, Tepsuper, Tmusuper] = ...
        solve3D_EigenEngine_Matrices(L0, wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml,s);
    toc
    %MultiCellReorder(A, b,xCells, yCells, N)
    %% parameters for a multicellreordering
    xCells = i;
    Grid = zeros(Nx, Ny, Nz);
    counter = 0;
    tic
    [permutedIndices, boundaryCells, interiorCells, intStor, hpart, vpart] = ...
    IndexPermutation3DOneDChain(N,xCells);
    toc
    %% visualize partitioning
    figure;
    plot3(boundaryCells(:,1), boundaryCells(:,2), boundaryCells(:,3), '.', 'markersize', 16);

    %% Execute a Symmetry Preserving Column Row Permutation Combination
    tic
    xind = zeros(3*Nx*Ny*Nz,1);
    yind = zeros(3*Nx*Ny*Nz,1);
    vals = ones(3*Nx*Ny*Nz,1);
    for j = 1:3*N(1)*N(2)*N(3)

       indexshift = permutedIndices(j);
       xind(j) = j;
       yind(j) = indexshift;
       %Q(i,indexshift) = 1;
    end
    Q = sparse(xind,yind,vals);
    toc
    %%test permutation matrix should permute index order to permuted indices

    %% Transform Equations
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    SymB = Q*b;

   %% Execute IterativeFourBlockSchur
   disp('four block schur')
   tic
   [Aschur, bmod, App, Apv, Avp, Avv, bp, bv] = ...
    FourblockSchur(SymA,SymB, hpart,vpart);
   toc

    %% solve equations
   iterations = 50:10: 300;
tol = 1e-14;
residuals = []; iterate0_array = []; iterateSchur_array = [];

[xUnred,flag1,relres1,iter1] = qmr(SymA,SymB,tol, length(SymB));% = A\b;
[xSchur,flag2,relres2,iter2] = qmr(Aschur,bmod,tol, length(bmod));

%% Extract Boundary from the unreduced System
xUnredBound = xUnred(1:hpart);
%% Construct the interior Solution from Schur Solution

recSolInt = Avv\(bv - Avp*xSchur);
tolSol = [xSchur; recSolInt];

    
end