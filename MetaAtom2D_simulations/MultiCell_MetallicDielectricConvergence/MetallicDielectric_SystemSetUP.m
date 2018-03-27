%% Multicell Single Partition Analysis with real dielectric + PML
close all
clear all

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 3.0;  % wavelength in L0
iterScaling = []; 
maxCellNumber = 1;
SingleCellSize = 80;
epsilon = -9.8+0.31i;
Sx = SingleCellSize; Sy = SingleCellSize;

%% ==================== Create Data Directory

for npml = [15, 0]
    k = 8;
    N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
    N0 = N; %record the original N for comparison
    Npml = [npml, npml];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
    xrange = k*[-1 1];  % x boundaries in L0
    yrange = k*[-1 1];  % y boundaries in L0

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
    featureDims = SingleCellSize*[1/4, 1/4, 1/6];
    [eps_air, cellIndices] =... 
        multiRandomCellDielectricSingleLayerSep(k, k, SingleCellSize, ...
        SingleCellSize, Npml,epsilon, featureDims); %% ADD coe to account for PML


    %% Set up the magnetic current source density.
    Mz = zeros(N);
    ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
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

    disp('reordering')
    tic
    [SymA, SymB, Q, permutedIndices, boundaryCells, interiorCells,pmlCell, ...
    hpart, vpart, pmlxpart, pmlypart, pmlCellDict] = ... 
        MultiCellBoundaryInteriorSLS_PML_partitioned(A, b, divx, divy, SingleCellSize, N,Npml);
    toc
    %% Perform an iterative schur complement: Aschur = sum(App - X^-1CX

    disp('schur complement')
    tic
    [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvcell, InvAvvStorage] = ...
        PrecompIterativeFourBlockSchurSingleSep_PartitionedPML(SymA, SymB, hpart, vpart,divx, divy,N,...
        SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart, pmlCellDict);
    toc
    x1 = A\b;
    x2 = Aschur\bmod;
    tolSol = ...
    MultiCellSchurInteriorSolGeneral(x2, Avvcell, ...
                            Avpcell, bvcell);     
    Hzrec = reshape(Q\tolSol, length(tolSol)^.5, length(tolSol)^.5);
    Hz = reshape(x1, Nx, Ny);
    fileName = strcat('Multicell2D_wvlen=',num2str(3),'Npml=',num2str(npml),...
        '_epsilon=',num2str(epsilon),'numCells=',num2str(k),'_SystemMatrices.mat');
    save(fileName,'Aschur', 'bmod', 'A', 'b', 'Q', 'Hz', 'Hzrec');
    
    
end


