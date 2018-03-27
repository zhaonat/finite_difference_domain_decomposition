%% Precomputed IterativeFourBlockSchur
% Four block as we have one eliminated domain (which consists of all of the
% subdomains in the metatoms, one reduced domain, and two coupling matrices
%% Parameter Descriptions

function [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell, bp, InvAvvStorage] = ...
    PrecompIterativeFourBlockSchur(SymA,SymB, hpart,vpart, xCells, yCells, N, cellIndices)

    %% cellIndices lists numbers from 1 to n indentifying what object is in each cell
    % cellindices is an xCellsxyCells array
    Nx = N(1); Ny = N(2);
    CDx = Nx/xCells; CDy = Ny/yCells;
    totCells = xCells*yCells;
    cellSeen = zeros(3,1);
    InvAvvStorage = cell(3,2);
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
    Aschur = App; bmod = bp;
    
    for i = 1:xCells*yCells
        cellType = cellIndices(i);
        if(cellSeen(cellType) == 0)
            cellSeen(cellType) = cellSeen(cellType)+1;
            Avv = Avvcell{i,i};
          
            [L,U] = lu(Avv);
            % store LU factorizations of Avv, which will be used later
            % to reconstruct the interior solutions
            InvAvvStorage{cellType, 1}= L;
            InvAvvStorage{cellType, 2} = U;

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
end