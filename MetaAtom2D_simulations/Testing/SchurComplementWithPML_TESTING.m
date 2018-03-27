%% Single layer Schur complement with PML setup
close all
clear

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.55*L0;  % wavelength in L0
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 40;
epsilon = 10;
Sx = SingleCellSize; Sy = SingleCellSize;

k = 1; xCells =k; yCells = k;

%% ================ Domain Size Parameters ==========================
N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
N0 = N;
Npml = [10 10];  
xrange = k*L0*[-1 1];  % x boundaries in L0
yrange = k*L0*[-1 1];  % y boundaries in L0

%% Generate Parameters for Domain with PML
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  

%% ==============Setup Dielectric ===================================
[eps_air, cellIndices] = ... 
    multiRandomCellDielectricSingleLayerSep(k, k, SingleCellSize,...
    SingleCellSize, Npml,epsilon); %% ADD coe to account for PML

%% setup point source 
Mz = zeros(N);
ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

Mz1 = zeros(N0);
Mz1(ceil(N0/2)) = 1;
%% ============ ==== FDFD Matrix Setup ==== ===========================%
%with pml
[A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices_unscaled(wvlen, xrange, ...
    yrange, eps_air, Mz, Npml);
    
%% Reorder the Matrix
[SymA, SymB, Q, permutedIndices, boundaryInd, interiorInd, pmlInd, hpart, vpart,...
    pmlxpart, pmlypart] = ... 
        MultiCellBoundaryInteriorSLS_PML(A, b, xCells, yCells, SingleCellSize, N,Npml);

%% Test of the PML Schur complement

    %% cellIndices lists numbers from 1 to n indentifying what object is in each cell
    %cellindices is an xCells x yCells array, so we may have to flatten
    %it
    cellIndices = cellIndices(:);
    %% add one extra element in cellindices to account for the pml domain
    cellIndices = [cellIndices; 4];
    Nx = N(1); Ny = N(2);
    CDx =SingleCellSize; CDy = SingleCellSize;
    totCells = xCells*yCells;
    cellSeen = zeros(4,1);
    InvAvvStorage = cell(3,2); %stores L and U in LU factorization
    %uppermost block = boundary cell partition
    % lower blocks = interior cell partition
    %% Partitioning of the SymA matrix into interface and non interface nodes
    Abound = SymA(1:hpart, 1:vpart);
    
    Aint = SymA(hpart+1:pmlxpart, vpart+1:pmlypart);
    Aib = SymA(hpart+1:pmlxpart,1:vpart);
    Abi = SymA(1:hpart, vpart+1:pmlypart);
    
    bBound = SymB(1:hpart);
    bint = SymB(hpart+1:pmlypart);
    bpml = SymB(pmlypart+1:end);
   
    %% PML SubBLOCK
    Apml = SymA(pmlxpart+1:end, pmlypart+1: end);
    Aci =  SymA(1:hpart, pmlxpart+1:end)
    Aic =  SymA(pmlxpart+1:end, 1:vpart)
    
    App = Abound;
    Avvcell = mat2cell(Aint, (CDx)^2*ones(totCells,1), (CDy)^2*ones(totCells,1));
    Apvcell = mat2cell(Abi, hpart, (CDx)^2*ones(totCells,1));
    Avpcell = mat2cell(Aib, (CDy)^2*ones(totCells,1), vpart);
    bvCell = mat2cell(bint, (CDx)^2*ones(totCells,1), 1);
    
    %% append the pml cells on
    bvCell{length(bvCell)+1} = bpml;
    Avvcell{totCells+1, totCells+1} = Apml;
    Apvcell{1,length(Apvcell)+1} = Aci;
    Avpcell{length(Avpcell)+1,1} = Aic;
    
    bp = bBound;
    Aschur = App;
    bmod = bp;
    
    %iterate through  thAve metaatoms...linear loop (x*y), but still N^2
    % loop is linear because the Avv's we want are all on a diagonal
    for i = 1:xCells*yCells+1 %add 1 because we have the pml domain
        cellType = cellIndices(i);
        if(cellSeen(cellType) == 0)
            cellSeen(cellType) = cellSeen(cellType)+1;
            Avv = Avvcell{i,i};
            [L,U] = lu(Avv);
            InvAvvStorage{cellType, 1}= L;
            InvAvvStorage{cellType, 2}= U;
        else
            cellSeen(cellType) = cellSeen(cellType)+1;
        end
        invAvvL = InvAvvStorage{cellType,1};
        invAvvU = InvAvvStorage{cellType,2};
        App = Abound;
        Apv = Apvcell{1, i};
        Avp = Avpcell{i, 1};
        bp = bBound;
        bv =bvCell{i};
        Comp = Apv*(invAvvU\(invAvvL\Avp));
        Aschur = Aschur - Comp; 
        bmod = bmod - Apv*(invAvvU\(invAvvL\bv));
       
    end


%% ============= Test Solve ==============================
x = A\b;
x2 = SymA\SymB;
boundaryEx = x2(1:hpart);
xschur = Aschur\bmod;
Hz = reshape(Q\x2, length(x)^.5, length(x)^.5);
imagesc(abs(Hz))


