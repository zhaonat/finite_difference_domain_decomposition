close all
clear;

%% Parameter Set up
%% Parameter Set-up
L0 = 1e-6;  % length unit: microns
wvlen = 5;  % wavelength in L0
for diel = [1];
cellsize = 20;
inclusionSize = cellsize-6;

for numcells = 1:1
    cellx = numcells; celly = numcells; cellz = 1;
    core = strcat('inclusionSize=', num2str(inclusionSize),'_numCells=', num2str(cellx*celly*cellz), +'_wvlen=', num2str(wvlen),'_diel=',num2str(diel),...
        'cellSize=',num2str(cellsize),'_3DMultiCell')
    Npml = [0 0 0];  % [Nx_pml Ny_pml]
    xrange = numcells*[-1 1];  % x boundaries in L0
    yrange = numcells*[-1 1];  % y boundaries in L0
    zrange = numcells*[-1 1];
    %[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  
    % domain is expanded to include PML

    %% Set up the permittivity.
    SingleCellDims = cellsize*[1, 1, 1];
    [eps_r, interiorCoords, borderCoords] = cubeDielectricGrid(cellx, celly,...
    cellz, SingleCellDims,...
    inclusionSize, diel);
    N = size(eps_r);

    %% Set up the current source density.
    Mz = zeros(N); My = Mz; Mx = Mz;
    ind_src = 5*[1 1 1];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
    JCurrentVector = [Mx; My; Mz];
    
    %% implement an iterative solve for all of the Schur complement elements
    %% Wonsoek's scalar parameter 1, -1, or 0
    s = -1;
    %=========================================
    tic
    [A,b,Ao, bo, omega, c0, Tepsuper, Tmusuper] = ...
        solve3D_EigenEngine_Matrices_dirichlet(L0, wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml,s);
    toc
    %MultiCellReorder(A, b,xCells, yCells, N)
    %% parameters for a multicellreordering
    tic
    [Q, permutedIndices, hpart, vpart, partitions] = ...
    IndexPermutationCubicLattice_SLS(N,borderCoords, interiorCoords)
    toc
  


    %% Transform Equations
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    SymB = Q*b;

   %% Execute IterativeFourBlockSchur
   disp('four block schur')
%    tic
%    [Aschur, bmod, App, Apv, Avp, Avv, bp, bv] = ...
%     FourblockSchur(SymA,SymB, hpart,vpart);
%    toc
   tic
   [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvcell] = ...
    ParallelIterativeFourBlockSchur3D(SymA,SymB, partitions);
    toc
    %% evaluating the factorization is on average, 100 times slower
    
end
end