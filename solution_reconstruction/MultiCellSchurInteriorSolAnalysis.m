%% Multicell interior reconstruction Analysis
close all

L0 = 1e-6;  % length unit: microns
wvlen = 5*L0;  % wavelength in L0
iterSchur = []; iterUnRed = [];

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
k = 2; %% Number of cells
    N = k*[SingleCellSize SingleCellSize];  % [Nx Ny]
    Npml = k*[0 0];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
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
    eps_air = multiRandomCellDielectric(k, k, Nx/k, Ny/k, Npml); %% ADD coe to account for PML
    
%     figure;
%     imagesc(abs(eps_air));

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

    [SymA, SymB, Q, permutedIndices, boundaryInd, interiorInd, hpart, vpart] = ... 
        MultiCellReorder(A, b, divx, divy, N);
    
    %% Perform an iterative schur complement: Aschur = sum(App - X^-1CX
    [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvcell] = ...
        IterativeFourBlockSchur(SymA, SymB, hpart, vpart,divx, divy,N);
    conditionNumbers = [conditionNumbers condest(Aschur) condest(A)];
    
    %% Now we need to establish formalism to extract the interior solutions
    %% for all cells
    xSchur = Aschur\bmod;
       
    %% Solver Code starts here
    xCells = divx; yCells = divy;
    %% break xSchur into cells
      
    tolSol = [xSchur];
    recSolIntArray = cell(xCells, yCells);
    %% Need a way to track all the correct cells to get the interior sol
    count = 1;
    for i = 1:xCells
        for j = 1:yCells
            Avv = Avvcell{count,count};
            bv = bvcell{count,1};
            Avp = Avpcell{count,1};
            yVec = bv - Avp*xSchur;
            recSolInt = Avv\yVec;
            tolSol = [tolSol; recSolInt];
            %this explicitly assuemes that all domains we Schur
            %complemented are square in shape...
            reshaped = reshape(recSolInt, SingleCellSize-2, SingleCellSize-2);
            figure;
            imagesc(abs(reshaped));
            recSolIntArray{i,j} = reshaped;
            count = count+1;
        end
    end
 figure
 imagesc(abs(reshape(A\b,Nx, Ny)));
    