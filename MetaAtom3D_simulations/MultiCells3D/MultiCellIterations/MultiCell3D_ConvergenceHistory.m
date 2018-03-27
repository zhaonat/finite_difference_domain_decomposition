close all
clear;

%% Parameter Set up
%% Parameter Set-up
L0 = 1e-6;  % length unit: microns
wvlen = 5;  % wavelength in L0
for diel = [-9.8+0.31i];
cellsize = 18;
inclusionSize = cellsize-6;

    for numcells = 3:3
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
            solve3D_EigenEngine_Matrices(L0, wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml,s);
        toc
        %MultiCellReorder(A, b,xCells, yCells, N)
        %% parameters for a multicellreordering
        tic
        [Q, permutedIndices, hpart, vpart, partitions] = ...
        IndexPermutationCubicLattice_SLS(N,borderCoords, interiorCoords);
        toc

    %     figure; 
    %     c2 = cell2mat(interiorCoords);
    %     scatter3(c2(:,1), c2(:,2), c2(:, 3), 'filled')
    %     hold on;
    %     scatter3(c2(:,4), c2(:,5), c2(:, 6), 'filled')
    % 
    %     hold on
    %     c3 = cell2mat(borderCoords);
    %     scatter3(c3(:,1), c3(:,2), c3(:, 3), 'filled')


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
        ParallelFourBlockSchur3D(SymA,SymB, partitions);
        toc

        %% solve equations
        tol = 1e-10;
        %% QMR
        disp('full solve');
        tic
        [xTol, flag2, relres2, iter2, resvec2] = ...
            qmr(A,b, 1e-10, length(b));
        toc
        
        disp('reduced solve')
        tic 
        [xSchur, flag1, relres1, iter1, resvec1] = ...
            qmr(Aschur, bmod, 1e-10, length(bmod));
        toc


        %% now we recover the residuals on the iterates

        save(strcat(core, '.mat'), 'b', 'bmod', 'resvec1', 'resvec2', ...
            'flag1', 'flag2', 'iter1', 'iter2', 'xSchur', 'xTol');


    end
end