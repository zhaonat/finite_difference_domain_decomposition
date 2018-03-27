%% Multicell Single Partition Analysis with real dielectric + PML
clear all

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1e8;  % wavelength in L0
maxCellNumber = 1;
SingleCellSize = 40;
epsilon = 1;
Sx = SingleCellSize; Sy = SingleCellSize;

%% ==================== Create Data Directory
condData = []; iterData = [];
for epsilon = [1,12]
iterScaling = []; conditions = [];

for numCells = 1:5
    N = [numCells*SingleCellSize+numCells+1 numCells*SingleCellSize+numCells+1];  % [Nx Ny]
    N0 = N; %record the original N for comparison
    Npml = 0*[5 5];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
    xrange = numCells*[-0.5  0.5];  % x boundaries in L0
    yrange = numCells*[-0.5, 0.5];  % y boundaries in L0

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
    FeatureDims = [SingleCellSize/4, SingleCellSize/4, SingleCellSize/6];
    [eps_air, cellIndices] =... 
        multiRandomCellDielectricSingleLayerSep(numCells, numCells, SingleCellSize, SingleCellSize, Npml,epsilon,FeatureDims); %% ADD coe to account for PML


    %% Set up the magnetic current source density.
    Mz = zeros(N);
    ind_src = [ceil(N/5) ceil(N/5)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;
    %Mz(75, 75) = 1;
    %scale = 1e-13; %% make matrix scaling look nicer
    [A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices(L0,wvlen, xrange, yrange, eps_air, Mz, Npml);
    %A = scale*A; b = scale*b;
    
    %% Set up the MultiCell Regime
    divx = numCells; divy = numCells;
    xCells = numCells; yCells = numCells; %% these specify the number of cells we will be dividing into
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

    conditions = [conditions; eigs(Aschur, 1, 'sm'), eigs(A, 1, 'sm')];

    
end
condData = [condData, conditions];
%iterData = [iterData, iterScaling];
figure; 
plot(conditions);
    
end
figure; 
plot(condData)