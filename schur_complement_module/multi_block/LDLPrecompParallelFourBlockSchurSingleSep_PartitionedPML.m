function [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell, InvAvvStorage] = ...
    PrecompParallelFourBlockSchurSingleSep_PartitionedPML(SymA, SymB, hpart, vpart, ...
    xCells, yCells, N, SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart, pmlCellDict)

    %% cellIndices lists numbers from 1 to n indentifying what object is in each cell
    %cellindices is an xCells x yCells array, so we may have to flatten
    %it
    cellIndices = cellIndices(:);
    
    Nx = N(1); Ny = N(2);
    CDx =SingleCellSize; CDy = SingleCellSize;
    totCells = xCells*yCells;
    InvAvvStorage = cell(3,2); %stores L and U in LU factorization
    numCells = totCells + (xCells+yCells)+1;

    
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
    Aci =  SymA(1:hpart, pmlypart+1:end);
    Aic =  SymA(pmlxpart+1:end, 1:vpart);
    
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
    SchurComplements = cell(numCells);
    SchurSources = cell(numCells);
    %iterate through  thAve metaatoms...linear loop (x*y), but still N^2
    % loop is linear because the Avv's we want are all on a diagonal
    for i = 1:numCells %add 1 because we have the pml domain
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
    end
    %% InvAVvstorage is sent to every one of the workers...but it will necessarily
    %% be small, it has only 3 LU factorizations
    parfor i = 1:numCells
        cellType = cellIndices(i);
        invAvvL = InvAvvStorage{cellType,1};
        invAvvU = InvAvvStorage{cellType,2};
        Apv = Apvcell{1, i};
        Avp = Avpcell{i, 1};
        bv =bvCell{i,1};
        Comp = Apv*(invAvvU\(invAvvL\Avp));
        SchurComplements{i} = Comp;
        SchurSources{i} = Apv*(invAvvU\(invAvvL\bv));
    end
    Aschur = App; bmod = bp;
    %% Construct final Schur complement
    for i = 1:numCells
       Aschur = Aschur - SchurComplements{i};
       bmod = bmod - SchurSources{i};
    end
    
end