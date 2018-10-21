%% Multicell interior reconstruction Analysis
close all
clear all

L0 = 1e-6;  % length unit: microns
wvlen = 5*L0;  % wavelength in L0
iterSchur = []; iterUnRed = [];
dielConst = 12;
SingleCellSize = 90;
Sx = SingleCellSize; Sy = SingleCellSize;


%% Create Data Directory
% cd(strcat('D:\Nathan\Documents\StanfordYearOne\Fan Group\FDFDSchurComplement\MetaAtom2D\MultiCells\'));
% 
% dataDir = strcat('MultiCell2D with SingleCellSize = ',int2str(SingleCellSize),...
%     'wvlen=',num2str(wvlen)));
% mkdir(dataDir);
% cd(strcat('D:\Nathan\Documents\StanfordYearOne\Fan Group\FDFDSchurComplement\MetaAtom2D\MultiCells\',dataDir));
% 
%     
conditionNumbers = [];
k = 1; %% Number of cells
xCells = k; yCells = k;
    N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
    Npml = k*[5 5];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
    xrange = k*L0*[-5 5];  % x boundaries in L0
    yrange = k*L0*[-5 5];  % y boundaries in L0

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
    [eps_air, cellIndices] = ...
    multiRandomCellDielectricSingleLayerSep(xCells, yCells, SingleCellSize,...
        SingleCellSize, Npml, dielConst)

    %% Set up the magnetic current source density.
    Mz = zeros(N);
    ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;
    %Mz(75, 75) = 1;
    
    [Hz, Ex, Ey, A, omega,b, Sxf, Dxf,Dyf] = solveTE(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);

    %% Set up the MultiCell Regime
    divx = k; divy = k; %% these specify the number of cells we will be dividing into
    CellDimx = N(1)/divx;
    CellDimy = N(2)/divy;
    
    %% Partitions
    [SymA, SymB, Q, permutedIndices, boundaryCells, interiorCells,pmlCell, ...
    hpart, vpart, pmlxpart, pmlypart] = ...
    MultiCellBoundaryInteriorSLS_PML(A, b, xCells, yCells,SingleCellSize, N, Npml);
    
    %% Perform an iterative schur complement: Aschur = sum(App - X^-1CX
    [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvcell, InvAvvStorage] = ...
    PrecompIterativeFourBlockSchurSingleSep_PML(SymA, SymB, hpart, vpart, ...
    xCells, yCells, N, SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart);
    
    %% Now we need to establish formalism to extract the interior solutions
    %% for all cells
    xSchur = Aschur\bmod;
       
    %% Solver Code starts here
    xCells = divx; yCells = divy;
    %% break xSchur into cells
      
    tolSol = [xSchur];
    %% Need a way to track all the correct cells to get the interior sol
    d = size(InvAvvStorage);
    for count = 1:length(bvcell)
        Avv = Avvcell{count,count};
        bv = bvcell{count};
        Avp = Avpcell{count,1};
        yVec = bv - Avp*xSchur;
        recSolInt = Avv\yVec;
        tolSol = [tolSol; recSolInt];
        %this explicitly assuemes that all domains we Schur
        %complemented are square in shape...

    end
    
 figure
 imagesc(abs(reshape(A\b,Nx, Ny)));
 figure()
 plot(abs(tolSol));
 
    