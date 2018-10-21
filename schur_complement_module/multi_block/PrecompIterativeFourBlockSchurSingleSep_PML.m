%% Precomputed IterativeFourBlockSchur
function [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell, InvAvvStorage] = ...
    PrecompIterativeFourBlockSchurSingleSep_PML(SymA, SymB, hpart, vpart, ...
    xCells, yCells, N, SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart)
    

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
    Aci =  SymA(1:hpart, pmlxpart+1:end);
    Aic =  SymA(pmlxpart+1:end, 1:vpart);
    
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
            tic
            [L,U] = lu(Avv);
            toc
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

end