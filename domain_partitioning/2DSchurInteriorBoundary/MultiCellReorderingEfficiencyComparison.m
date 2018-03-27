clear
close all
%% Note Schur PML currently does not work in the single cell case
%% === ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 90;
epsilon = 1;
Sx = SingleCellSize; Sy = SingleCellSize;
   

%% REMEMBER WE CANNOT ANALYZE LARGE MATRICES WITH THIS
maxEigens = []; conditionNums = [];
unreducedEigenSpectra = [];
reducedEigenSpectra = [];

for k = 6:6
   
    wvlen = 2.0;  % wavelength in L0

    N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
    Npml = k*[0 0];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
%     xrange = k*L0*[-1 1]*.5;  % x boundaries in L0
%     yrange = k*L0*[-1 1]*.5;  % y boundaries in L0
    xrange = k*[-1 1]*0.5;  % x boundaries in L0
    yrange = k*[-1 1]*0.5;  % y boundaries in L0

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
    ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;
    %Mz(75, 75) = 1;
    %scale = 1e-20; %% make matrix scaling look nicer
    [A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
    %A = scale*A; b = scale*b;


    %% Set up the MultiCell Regime
    divx = k; divy = k;
    xCells = k; yCells = k; %% these specify the number of cells we will be dividing into
    CellDimx = N(1)/divx;
    CellDimy = N(2)/divy;
    disp('reordering')
    tic
    [SymA, SymB, Q, permutedIndices, boundaryCells, interiorCells,pmlCell, ...
    hpart, vpart, pmlxpart, pmlypart] = ... 
        MultiCellBoundaryInteriorSLS_PML(A, b, divx, divy, ...
        SingleCellSize, N, Npml);
    toc
    
    %% now try the vectorized interior exterior version
    bounds = [200, 200, 400, 400];
    tic
    [SymA, SymB, Q, permutedIndices, hpart, vpart, map] =...
    RectangularReorder(A,b, N, bounds) ;
    toc
    %% this is obviously the fastest reordering, as it is vectorized
end