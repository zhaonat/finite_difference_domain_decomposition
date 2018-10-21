%% Multicell Single Partition Analysis with real dielectric + PML
close all
clear 

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.55*L0;  % wavelength in L0
maxCellNumber = 1;
SingleCellSize = 80;
epsilon = 12;
Sx = SingleCellSize; Sy = SingleCellSize;


rediterScaling = []; 
unrediterScaling = [];
cellSizes = 1:2
dims = [];
for k = 1:10
    g = k;
    dims = [dims, g*k];
    disp(dims)
    N = [g*SingleCellSize+g+1 k*SingleCellSize+k+1];  % [Nx Ny]
    N0 = N; %record the original N for comparison
    Npml = [25 25];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
    xrange = g*L0*[-1 1];  % x boundaries in L0
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
    featureDims = [SingleCellSize/2, SingleCellSize/4, SingleCellSize/6];
    [eps_air, cellIndices] =... 
        multiRandomCellDielectricSingleLayerSep(k, g, ...
        SingleCellSize, SingleCellSize, Npml,epsilon, featureDims); %% ADD coe to account for PML


    %% Set up the magnetic current source density.
    Mz = zeros(N);
    ind_src = [Npml(1)+10 Npml(2)+10];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;
    %Mz(75, 75) = 1;
    %scale = 1e-13; %% make matrix scaling look nicer

    [A, omega,b, Sxf, Syf, Sxb, Syb] = solveTE_Matrices_unscaled(wvlen, xrange, yrange, eps_air, Mz, Npml);
    %A = scale*A; b = scale*b;

    [Pl, Pr] = SCSymmetrizer2D(Sxf, Syf, Sxb, Syb);
    A = Pl^-1*(A*Pr^-1);
    b = Pl\b;
    
    %% Set up the MultiCell Regime
    divx = g; divy = k;
    xCells = g; yCells = k; %% these specify the number of cells we will be dividing into
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

    %% Now we need to establish formalism to extract the interior solutions

    %% perform iterations Analysis
    tol = 1e-14;
    maxitun = length(b);
    maxitred = length(bmod);
    %% We need to extract the convergence history for the border and the interior
    disp('full solve')
    tic 
    [x1, flag1, relres1, iter1, resvec1] = qmr(SymA, SymB, tol, maxitun);
    toc
    flag1
    iter1
    disp('reduced solve')
    tic
    [x2, flag2, relres2, iter2, resvec2] = qmr(Aschur, bmod, tol, maxitred);
    toc
    flag2
    iter2
    rediterScaling = [rediterScaling; iter2]
    %the unreduced case is problematic as their is significant stagnation
    if(flag1 == 1)
        unrediterScaling = [unrediterScaling; length(b)];
    else
        unrediterScaling = [unrediterScaling; iter1];
    end
end

%% Plot cell scaling properties
g = figure;
plot(dims, rediterScaling);
hold on
plot(dims, unrediterScaling);
title('Cell Scaling for PartitionedPML (fully complemented out)')
xlabel('Iteration #')
ylabel('log_{10}(Relative R1esidual)')
saveas(g, strcat('Npml=',num2str(Npml(1)), ...
    'IterationScalingWithPartitionedPML (fully complemented out).fig'))




