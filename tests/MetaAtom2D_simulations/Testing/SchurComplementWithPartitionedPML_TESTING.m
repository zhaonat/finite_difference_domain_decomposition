%% Single layer Schur complement with Partitioned PML setup
close all
clear

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.55*L0;  % wavelength in L0
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 50;
epsilon = 12;
Sx = SingleCellSize; Sy = SingleCellSize;

k = 3; xCells =k; yCells = k;

%% ================ Domain Size Parameters ==========================
N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
N0 = N;
Npml = [5 5];  
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

%% ============ ==== FDFD Matrix Setup ==== ===========================%
%with pml
[A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices_unscaled(wvlen, xrange, ...
    yrange, eps_air, Mz, Npml);
    
%% Reorder the Matrix
[SymA, SymB, Q, permutedIndices, boundaryInd, interiorInd, pmlInd, hpart, vpart,...
    pmlxpart, pmlypart, pmlCellDict] = ... 
        MultiCellBoundaryInteriorSLS_PML_partitioned(A, b, xCells, yCells, SingleCellSize, N,Npml);

%% Test of the PML Schur complement

    %% cellIndices lists numbers from 1 to n indentifying what object is in each cell
    %cellindices is an xCells x yCells array, so we may have to flatten
    %it
    cellIndices = cellIndices(:);
    
    Nx = N(1); Ny = N(2);
    CDx =SingleCellSize; CDy = SingleCellSize;
    totCells = xCells*yCells;
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
    Aci =  SymA(1:hpart, pmlypart+1:end)
    Aic =  SymA(pmlxpart+1:end, 1:vpart)
    
    App = Abound;
    Avvcell = mat2cell(Aint, (CDx)^2*ones(totCells,1), (CDy)^2*ones(totCells,1));
    Apvcell = mat2cell(Abi, hpart, (CDx)^2*ones(totCells,1));
    Avpcell = mat2cell(Aib, (CDy)^2*ones(totCells,1), vpart);
    bvCell = mat2cell(bint, (CDx)^2*ones(totCells,1), 1);
    
    
    %% append the pml cells on (this time, it's a set of partitioned pml parts)
    pmllargeCell = size(pmlCellDict{1}); % the edge slicess
    pmlsmallCell = size(pmlCellDict{end}); %the corners
    perimeter = 2*xCells+2*yCells;
    
    s1 = 2*pmllargeCell(1)*ones(perimeter/2, 1); %partition for larger pmls
    s2 = 4*pmlsmallCell(1); 
    sizes = [s1;s2];
    Apmlcell = mat2cell(Apml, sizes, sizes);
    bpmlcell = mat2cell(bpml, sizes, 1);
    Avppml = mat2cell(Aic, sizes, hpart);
    Apvpml = mat2cell(Aci, vpart, sizes);
    Apvcell = horzcat(Apvcell, Apvpml);
    Avpcell = vertcat(Avpcell, Avppml);
    ind = length(Avvcell)
    for i = 1:length(Apmlcell)
        Avvcell{ind+i, ind+i} = Apmlcell{i,i}; 
        %ALSO modify cellindices to account for the partitioned pml domain
        cellIndices = [cellIndices; 3+i];%start with 3 as 1-3 specify interior dielectric shapes
    end
    %% seen dictionary over here because this is where cellIndices is at final length
    cellSeen = zeros(max(cellIndices),1);

    
    %% Now concatenate the cell arrays of the partitioned pml onto the overall
    %% cell array
    bvCell = vertcat(bvCell,bpmlcell);
    bp = bBound;
    Aschur = App;
    bmod = bp;
    
    %iterate through  thAve metaatoms...linear loop (x*y), but still N^2
    % loop is linear because the Avv's we want are all on a diagonal
    for i = 1:xCells*yCells+xCells +yCells +1 %add 1 because we have the pml domain
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
        bv =bvCell{i,1};
        Comp = Apv*(invAvvU\(invAvvL\Avp));
        Aschur = Aschur - Comp; 
        bmod = bmod - Apv*(invAvvU\(invAvvL\bv));
       
    end

%% Visualize All diagonal cells
% d2 = size(Avvcell);
% nnzsum = 0;
% for i = 1:d2[1]
%    figure()
%    spy(Avvcell{i,i});
%    nnzsum = nnzsum + nnz(Avvcell{i,i})
% end


%% ============= Test Solve ==============================
x = A\b;
x2 = SymA\SymB;
boundaryEx = x2(1:hpart);
xschur = Aschur\bmod;
Hz = reshape(Q\x2, length(x)^.5, length(x)^.5);
tolSol =MultiCellSchurInteriorSolGeneral(xschur, Avvcell, InvAvvStorage, ...
    Avpcell, bvCell);
Hzrec = reshape(Q\tolSol, length(tolSol)^.5, length(tolSol)^.5);
imagesc(abs(Hz))

figure()
plot(abs(x2))
hold on
plot(abs(xschur))

