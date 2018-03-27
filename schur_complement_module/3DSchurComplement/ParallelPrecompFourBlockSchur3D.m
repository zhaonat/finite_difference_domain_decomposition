%% Precomputed IterativeFourBlockSchur
function [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell] = ...
    ParallelPrecompFourBlockSchur3D(SymA,SymB, partitions)
    
    %% cellIndices lists numbers from 1 to n indentifying what object is in each cell
        %cellindices is an xCellsxyCells array
    totCells = length(partitions(2:end))
    subsystemSize = partitions(2);
    hpart = partitions(1); vpart = hpart;
    %uppermost block = boundary cell partition
    % lower blocks = interior cell partition
    Abound = SymA(1:hpart, 1:hpart);
    

    Aint = SymA(hpart+1:end, vpart+1:end);
    Aib = SymA(hpart+1:end,1:vpart);
    Abi = SymA(1:hpart, vpart+1:end);
    
    bBound = SymB(1:hpart);
    bint = SymB(hpart+1:end);
   
    App = Abound;
    Avvcell = mat2cell(Aint, subsystemSize*ones(totCells,1), subsystemSize*ones(totCells,1));

    Apvcell = mat2cell(Abi, hpart, subsystemSize*ones(totCells,1));
    Avpcell = mat2cell(Aib, subsystemSize*ones(totCells,1), vpart);
    bvCell = mat2cell(bint, subsystemSize*ones(totCells,1), 1);
    bp = bBound;
    Aschur = App;
    bmod = bp;
     
    %iterate through  thAve metaatoms...linear loop (x*y), but still N^2
    % loop is linear because the Avv's we want are all on a diagonal
    SchurComplements = cell(totCells+1); %% this stores all the BD^-1C
    SchurSources = cell(totCells+1);
    InvAvvL = cell(1);
    InvAvvD = cell(1);
    InvAvvP = cell(1);
    %% FOR NOW, ALL CELL INCLUSIONS ARE THE SAME
    [L,D,P] = ldl(Avvcell{1,1});
    
    parfor i = 1:totCells
        Avv = Avvcell{i,i};

        Apv = Apvcell{i};
        Avp = Avpcell{i};
        bv =bvCell{i};
        Comp = Apv*(P.'\(L.'\(D\(L\(P\Avp)))));
        SchurComplements{i} = Comp;
        %yv = qmr(Avv, bv, 1e-8, length(bv));
        SchurSources{i} = Apv*(P.'\(L.'\(D\(L\(P\bv)))));
    end
    %% Construct final Schur complement
    for i = 1:totCells
       Aschur = Aschur - SchurComplements{i};
       bmod = bmod-SchurSources{i};
    end
    
end