%% Precomputed IterativeFourBlockSchur
function [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell, InvAvvStorage] = ...
    ParallelIterativeFourBlockSchurSingleSep_PML(SymA, SymB, hpart, vpart, ...
    xCells, yCells, N, SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart)
    

    %% cellIndices lists numbers from 1 to n indentifying what object is in each cell
    %cellindices is an xCells x yCells array, so we may have to flatten
    %it
    cellIndices = cellIndices(:);
    %% add one extra element in cellindices to account for the pml domain
    cellIndices = [cellIndices; 4];
    CDx =SingleCellSize; CDy = SingleCellSize;
    totCells = xCells*yCells;
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
    SchurComplements = cell(xCells*yCells+1); %% this stores all the BD^-1C
    SchurSources = cell(xCells*yCells+1);
    Avvflat = cell(xCells*yCells+1);
    %very simple transformation of Avvcell to a flat array
    for i = 1:xCells*yCells+1 %+1 is for the PML
        Avvflat{i} = Avvcell{i,i};
    end
    InvAvvStorageL = cell(totCells+1);
    InvAvvStorageU = cell(totCells+1);
    parfor i = 1:xCells*yCells+1
        Avv = Avvflat{i};
        [L,U] = lu(Avv);
        InvAvvStorageL{i}= L;
        InvAvvStorageU{i}= U;
        Apv = Apvcell{1, i};
        Avp = Avpcell{i, 1};
        bv =bvCell{i};
        Comp = Apv*(U\(L\Avp));
        SchurComplements{i} = Comp;
        SchurSources{i} = Apv*(U\(L\bv));
    end
    InvAvvStorage = [InvAvvStorageL.', InvAvvStorageU.'];
    %% Construct final Schur complement
    for i = 1:xCells*yCells+1
       Aschur = Aschur - SchurComplements{i};
       bmod = bmod-SchurSources{i};
    end

end