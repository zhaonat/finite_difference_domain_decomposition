%% Multicell Single Partition Analysis with real dielectric + PML
close all
clear all

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.55*L0;  % wavelength in L0
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 90;
epsilon = 12i;
Sx = SingleCellSize; Sy = SingleCellSize;

%% ==================== Create Data Directory
directory = strcat('D:\Nathan\Documents\StanfordYearOne\Fan Group\FDFDSchurComplement\MetaAtom2D\',...
                    'MultiCellsSingleLayerSeparator_PML');

for k = 1:maxCellNumber
    N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
    N0 = N; %record the original N for comparison
    Npml = k*[10 10];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
    xrange = k*L0*[-1 1];  % x boundaries in L0
    yrange = k*L0*[-1 1];  % y boundaries in L0

    %% Note on grid resolution of the system
    % dx/(wvlen) ~1/20 or smaller

    [xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
    %FINAL GRID PARAMETERS ARE DETERMINED AT THIS POINT
    %xrange and yrange are slightly larger
    %% Cell Division Setup
    M= N(1)*N(2);

    %% NOTE dL is not in SI units when it comes out of domain_with_pml; for our purposes, that is okay
    %for now
    resolutionFactor = max([dL(1)/N(1) dL(2)/N(2)]); %dx/N ~meters
    %spatially, what is the smallestlength scale that we have to resolve
    Nx = N(1); Ny = N(2);

    %% Set up the permittivity.
    [eps_air, cellIndices] =... 
        UniformLayerCellDielectricSingleLayerSep(k, k, SingleCellSize, SingleCellSize, Npml,epsilon); %% ADD coe to account for PML


    %% Set up the magnetic current source density.
    Mz = zeros(N);
    ind_src = [ceil(N/5) ceil(N/5)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;
    %Mz(75, 75) = 1;
    %scale = 1e-13; %% make matrix scaling look nicer
    [A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices_unscaled(wvlen, xrange, yrange, eps_air, Mz, Npml);
    %A = scale*A; b = scale*b;
    
    %% Set up the MultiCell Regime
    divx = k; divy = k;
    xCells = k; yCells = k; %% these specify the number of cells we will be dividing into
    CellDimx = N(1)/divx;
    CellDimy = N(2)/divy
    
    disp('unpartitioned reordering')
    tic
    [SymA1, SymB1, Q, permutedIndices, boundaryCells, interiorCells,pmlCell, ...
    hpart, vpart, pmlxpart, pmlypart] = ... 
        MultiCellBoundaryInteriorSLS_PML(A, b, divx, divy, SingleCellSize, N,Npml);
    unpr = toc
    %% Perform an iterative schur complement: Aschur = sum(App - X^-1CX
    
    disp('unpartitioned schur complement')
    tic
    [Aschur1, bmod1, App, Avvcell, Apvcell, Avpcell, bvcell, InvAvvStorage] = ...
        PrecompIterativeFourBlockSchurSingleSep_PML(SymA1, SymB1, hpart, vpart,divx, divy,N,...
        SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart);
    unpsc = toc
    
    disp('partitioned reordering')
    tic
    [SymA, SymB, Q, permutedIndices, boundaryCells, interiorCells,pmlCell, ...
    hpart, vpart, pmlxpart, pmlypart, pmlCellDict] = ... 
        MultiCellBoundaryInteriorSLS_PML_partitioned(A, b, divx, divy, SingleCellSize, N,Npml);
    pr = toc
    %% Perform an iterative schur complement: Aschur = sum(App - X^-1CX
    
    disp('partitioned schur complement')
    tic
    [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvcell, InvAvvStorage] = ...
        PrecompIterativeFourBlockSchurSingleSep_PartitionedPML(SymA, SymB, hpart, vpart,divx, divy,N,...
        SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart, pmlCellDict);
    prsc = toc

    %% Aschur resvecs
    [x1, flag1, relres1, iter1, resvec1] = qmr(Aschur1, bmod1, 1e-10, length(b));
    [x2, flag2, relres2, iter2, resvec2] = qmr(Aschur, bmod, 1e-10, length(b));
    
    tolSol = ...
    MultiCellSchurInteriorSolGeneral(x2, Avvcell, InvAvvStorage, ...
                            Avpcell, bvcell);     
    Hzrec = reshape(Q\tolSol, length(tolSol)^.5, length(tolSol)^.5)
    visabs(Hzrec, xrange, yrange);
    figure()
    semilogy(resvec1)
    hold on
    semilogy(resvec2)
    legend('unpartitioned', 'partitioned')
    
    
end


