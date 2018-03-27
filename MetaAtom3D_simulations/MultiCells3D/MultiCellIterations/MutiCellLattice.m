close all
clear;

%% Parameter Set up
%% Parameter Set-up
L0 = 1e-6;  % length unit: microns
wvlen = 6 ;  % wavelength in L0
diel = 4i;

for i = 2:2
    cellx = i; celly = i; cellz = 2;
    core = strcat('numCells=', num2str(cellx*celly*cellz), +'_wvlen=', num2str(wvlen),'_diel=',num2str(diel),...
        '_3DMultiCell')
    Npml = [0 0 0];  % [Nx_pml Ny_pml]
    xrange = i*[-1 1];  % x boundaries in L0
    yrange = i*[-1 1];  % y boundaries in L0
  
    zrange = i*[-1 1];
    %[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  
    % domain is expanded to include PML
    %% Set up the permittivity.
    inclusionSize = 4;
    SingleCellDims = 7*[1,1,1];
    [eps_r, interiorCoords, borderCoords] = cubeDielectricGrid(cellx, celly,...
    cellz, SingleCellDims,...
    inclusionSize, diel);
    N = size(eps_r);

    %% Set up the current source density.
    Mz = zeros(N); My = Mz; Mx = Mz;
    ind_src = [5 5 5];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
    JCurrentVector = [Mx; My; Mz];

    %% Wonsoek's scalar parameter 1, -1, or 0
    s = 0;
    %=========================================
    tic
    [A,b,Ao, bo, omega, c0, Tepsuper, Tmusuper] = ...
        solve3D_EigenEngine_Matrices(L0, wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml,s);
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
   %% fast ldl factorization is inaccurate for our super small cell sizes
   [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell] = ...
    ParallelIterativeFourBlockSchur3D(SymA,SymB, partitions);
    toc
    
    %% solve equations
    tol = 1e-6;
    disp('unreduced')
    tic
    [xUnred,flag1,relres1,iter1, resvec1] = qmr(SymA,SymB,tol, length(SymB));% = A\b;
    toc
    disp('reduced')
    tic
    [xSchur,flag2,relres2,iter2, resvec2] = qmr(Aschur,bmod,tol, length(bmod));
    toc
    f = figure;
    semilogy(resvec1/norm(b))
    hold on;
    semilogy(resvec2/norm(bmod))
    legend('unreduced', 'reduced')
    saveas(f, strcat(core, '_Convhist.fig'));
    save(strcat(core, '.mat'), 'b', 'bmod', 'resvec1', 'resvec2', 'A', 'Aschur');
    
    
end