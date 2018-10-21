%% Multicell Single Partition Analysis
%close all
clear

%% === ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 3;  % wavelength in L0
iterSchur = []; iterUnRed = [];
maxCellNumber = 15  ;
SingleCellSize = 80;
epsilon = 12;
Sx = SingleCellSize; Sy = SingleCellSize;

Total = [];
for pml = [0,15]
conditionNums = [];
    for k = 16:24
        N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
        Npml = [pml, pml];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
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
        xCells = k; yCells = k;
        featureDims = [SingleCellSize/4, SingleCellSize/4, SingleCellSize/6];
        celldimx = SingleCellSize; celldimy = SingleCellSize;
           [eps_air, cellIndices] = ...
        multiRandomCellDielectricSingleLayerSep(xCells, yCells,celldimx, ...
        celldimy, Npml, epsilon, featureDims);
        %% Set up the magnetic current source density.
        Mz = zeros(N);
        ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
        Mz(ind_src(1), ind_src(2)) = 1;
        %Mz(75, 75) = 1;
        scale = 1e-13; %% make matrix scaling look nicer
        [A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
        A = scale*A; b = scale*b;

        %% Set up the MultiCell Regime
        divx = k; divy = k;
        xCells = k; yCells = k; %% these specify the number of cells we will be dividing into
        CellDimx = N(1)/divx;
        CellDimy = N(2)/divy;
        disp('reordering')
        tic
         [SymA, SymB, Q, permutedIndices, boundaryCells, interiorCells,pmlCell, ...
    hpart, vpart, pmlxpart, pmlypart,pmlCellDict, pmlBoundaryCellCoords] = ...
    MultiCellBoundaryInteriorSLS_PML_partitioned(A, b, xCells, yCells,SingleCellSize, N, Npml);
        toc
        %% Perform an iterative schur complement: Aschur = sum(App - X^-1CX

        disp('schur complement')
        tic
        [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell, InvAvvStorage, Aint] = ...
         PrecompIterativeFourBlockSchurSingleSep_PartitionedPML(SymA, SymB, hpart, vpart, ...
            xCells, yCells, N, SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart, pmlCellDict);
        toc
        %conditionNumbers = [conditionNumbers condest(Aschur) condest(A)];

        %% Compare condition Numbers
        tic
        conditionNums = [conditionNums;k,  condest(App), condest(A), condest(Aschur)]
        toc

    end
    Total = [Total; conditionNums];

    figure()
    semilogy(conditionNums(:,2:end), 'linewidth', 2);
    xlabel('Total Number of Cells')
    title(strcat('ConditionData_PML=',num2str(Npml(1))));
    ylabel('Condition Number')
    legend('App', 'A', 'Aschur')
    grid()
    condApp = conditionNums(:,2);
    condA = conditionNums(:,3);
    condSchur = conditionNums(:,4);
    cellSizes = conditionNums(:,1)
    save(strcat('BigConditionData_PML=',num2str(Npml(1)),'.mat'),...
        'condApp','condA', 'condSchur','cellSizes', 'epsilon', 'wvlen')

end