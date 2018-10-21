%% Precomputed Iterative Four Block Schur

%% Multicell TIMING ANALYSIS: dual LAYER Partition
close all
clear;

L0 = 1e-6;  % length unit: microns
wvlen = 1.550;  % wavelength in L0
iterSchur = []; iterUnRed = [];

SingleCellSize = 80;
Sx = SingleCellSize; Sy = SingleCellSize;


    
conditionNumbers = [];
 k = 3;
    N = k*[SingleCellSize SingleCellSize];  % [Nx Ny]
    xCells = k; yCells = k;
    Npml = k*[0 0];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
    xrange = k*[-2 2];  % x boundaries in L0
    yrange = k*[-2 2];  % y boundaries in L0

    %% Note on grid resolution of the system
    % dx/(wvlen) ~1/20 or smaller

    [xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
    % FINAL GRID PARAMETERS ARE DETERMINED AT THIS POINT
    % xrange and yrange are slightly larger
    %% Cell Division Setup
    M= N(1)*N(2);

    %% NOTE dL is not in SI units when it comes out of domain_with_pml; for our purposes, that is okay
    %for now
    resolutionFactor = max([dL(1)/N(1) dL(2)/N(2)]); %dx/N ~meters
    %spatially, what is the smallestlength scale that we have to resolve
    Nx = N(1); Ny = N(2);

    %% Set up the permittivity.
    [eps_air, cellIndices] = multiRandomCellDielectric(xCells, yCells,Nx/k, Ny/k, Npml); %% ADD coe to account for PML
    
    figure;
    imagesc(abs(eps_air));

    %% Set up the magnetic current source density.
    Mz = zeros(N);
    ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;
    %Mz(75, 75) = 1;
    
    [A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
%     A = 1e-19*A;
%     b = 1e-19*b;
    %% Set up the MultiCell Regime
    divx = k; divy = k; %% these specify the number of cells we will be dividing into
    CellDimx = N(1)/divx;
    CellDimy = N(2)/divy;
    
    disp('cell reorder')
    tic
    [SymA, SymB, Q, permutedIndices, interiorInd,boundaryInd, hpart, vpart] = ... 
        MultiCellReorder(A, b, divx, divy, N);
    toc
    %% Comparison IterativeFourBlockSchur
    disp('it schur benchmark') %% RESULT: this is SEVERAL TIMES FASTER TO DO
    tic
     [Aschur1, bmod1, App1, Avvcell1, Apvcell1, Avpcell1, bvcell1] = ...
        PrecompIterativeFourBlockSchur(SymA, SymB, hpart, vpart,divx, divy,N, cellIndices);
    toc
    %% SOLVER CODE BEGINS HERE
    
    %% cellIndices lists numbers from 1 to n indentifying what object is in each cell
    %cellindices is an xCellsxyCells array
    %% indexing and cell-arrays
    disp('indexing and cell-arrays')
    tic
    Nx = N(1); Ny = N(2);
    CDx = Nx/xCells; CDy = Ny/yCells;
    totCells = xCells*yCells;
    cellSeen = zeros(3,1);
    LUAvvStorage = cell(3,1);
    %uppermost block = boundary cell partition
    % lower blocks = interior cell partition
    Abound = SymA(1:hpart, 1:vpart);
    

    Aint = SymA(hpart+1:end, vpart+1:end);
    Aib = SymA(hpart+1:end,1:vpart);
    Abi = SymA(1:hpart, vpart+1:end);
    
    bBound = SymB(1:hpart);
    bint = SymB(hpart+1:end);
   
    App = Abound;
    Avvcell = mat2cell(Aint, (CDx-2)^2*ones(totCells,1), (CDy-2)^2*ones(totCells,1));
    Apvcell = mat2cell(Abi, hpart, (CDx-2)^2*ones(totCells,1));
    Avpcell = mat2cell(Aib, (CDy-2)^2*ones(totCells,1), vpart);
    bvCell = mat2cell(bint, (CDx-2)^2*ones(totCells,1), 1);
    bp = bBound;
    Aschur = App;
    bmod = bp;
    toc
    disp('schur complements')
    tic
    for i = 1:xCells*yCells
            cellType = cellIndices(i);
            if(cellSeen(cellType) == 0)
                cellSeen(cellType) = cellSeen(cellType)+1;
                Avv = Avvcell{i,i};
                
                [L,U] = lu(Avv);
                
                LUAvvStorage{cellType, 1}= L;
                LUAvvStorage{cellType, 2}= U;
            else
                cellSeen(cellType) = cellSeen(cellType)+1;
            end
            L = LUAvvStorage{cellType,1};
            U = LUAvvStorage{cellType,2};
            App = Abound;
            Apv = Apvcell{1, i};
            Avp = Avpcell{i, 1};
            bp = bBound;
            bv =bvCell{i,1};
            disp('computing cost')
            tic
            Comp = Apv*(U\(L\Avp));
            Aschur = Aschur - Comp; 
            bmod = bmod - Apv*(U\(L\bv));
            toc

     end
    toc
    