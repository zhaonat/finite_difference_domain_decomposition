%% THIS FUNCTION CAN BE USED FOR BOTH SINGLE AND DUAL LAYER SEPARATORS!
%% the interior solution really just needs InvAvvStorage

function [tolSol, recSolIntArray] = ...
    MultiCellSchurInteriorSolLDL(xSchur, Avvcell, InvAvvStorage, Avpcell, bvcell,...
    xCells, yCells, cellIndices)
    
    %% Description of the Input Arguments
    %  xSchur - reduced solution
    %  [Avv, Avp; Apv, App]; Avv is the block that we complement out
    %  InvAvvStorage - cell containing LU factorizatoins of matrices in Avv
    %  Avpcell - couplings between Avv and App
    %  bvcell - source vector (Ax = b)
    %  xCell - number of metatoms in x direction
    %  yCells - number of metatoms in y direction
    %  SingleCellSize NxN nodes for each metatom
    %  cellindices
    %% Block Matrix form
    % Avv Avp /We reduce into Avv
    % Apv App /App represents the domain we get rid of
    
    tolSol = [xSchur]; %tolSol is not reordered 
    recSolIntArray = cell(xCells, yCells);
    %% Need a way to track all the correct cells to get the interior sol
    count = 1;
    for i = 1:xCells
        for j = 1:yCells
            
            L = InvAvvStorage{count, 1};
            D = InvAvvStorage{count, 2};
            P = InvAvvStorage{count, 3};
            bv = bvcell{count,1};
            Avp = Avpcell{count,1};
            yVec = bv - Avp*xSchur;
            recSolInt = P.'\(L.'\(D\(L\(P\yVec))));
            tolSol = [tolSol; recSolInt];
            reshaped = reshape(recSolInt, length(recSolInt)^.5, length(recSolInt)^.5);
            recSolIntArray{i,j} = reshaped;
            count = count+1;
        end
    end
    
    
end