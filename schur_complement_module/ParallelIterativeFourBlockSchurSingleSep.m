%% Precomputed IterativeFourBlockSchur
function [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell,...
    SchurComplements, SchurSources] = ...
    ParallelIterativeFourBlockSchurSingleSep(SymA, SymB, hpart, vpart, ...
    xCells, yCells, SingleCellSize)
    
    %% cellIndices lists numbers from 1 to n indentifying what object is in each cell
        %cellindices is an xCells x yCells array, so we may have to flatten
        %it
        
    CDx =SingleCellSize; CDy = SingleCellSize;
    totCells = xCells*yCells;

    %uppermost block = boundary cell partition
    % lower blocks = interior cell partition
    Abound = SymA(1:hpart, 1:vpart);
    
    
    Aint = SymA(hpart+1:end, vpart+1:end);
    Aib = SymA(hpart+1:end,1:vpart);
    Abi = SymA(1:hpart, vpart+1:end);
    
    bBound = SymB(1:hpart);
    bint = SymB(hpart+1:end);
   
    App = Abound;
    Avvcell = mat2cell(Aint, (CDx)^2*ones(totCells,1), (CDy)^2*ones(totCells,1));
    Apvcell = mat2cell(Abi, hpart, (CDx)^2*ones(totCells,1));
    Avpcell = mat2cell(Aib, (CDy)^2*ones(totCells,1), vpart);
    bvCell = mat2cell(bint, (CDx)^2*ones(totCells,1), 1);
    bp = bBound;
    Aschur = App;
    bmod = bp;
    %iterate through  the metaatoms...linear loop (x*y), but still N^2
    % loop is linear because the Avv's we want are all on a diagonal
    SchurComplements = cell(xCells*yCells);
    SchurSources = cell(xCells*yCells);
    Avvflat = cell(xCells*yCells);
    %very simple transformation of Avvcell to a flat array
    for i = 1:xCells*yCells
        Avvflat{i} = Avvcell{i,i};
    end
    parfor i = 1:xCells*yCells
        Avv = Avvflat{i};
        [L,U] = lu(Avv);
        Apv = Apvcell{1, i};
        Avp = Avpcell{i, 1};
        bv =bvCell{i,1};
        Comp = Apv*(U\(L\Avp));
        SchurComplements{i} = Comp;
        SchurSources{i} = Apv*(U\(L\bv));
    end
    %% Construct final Schur complement
    for i = 1:xCells*yCells
       Aschur = Aschur - SchurComplements{i};
       bmod = bmod-SchurSources{i};
    end
end