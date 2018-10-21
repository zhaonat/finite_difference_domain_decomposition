%% THIS FUNCTION, IN THEORY,
%% CAN BE USED FOR BOTH SINGLE AND DUAL LAYER SEPARATORS

function [tolSol, recSolIntArray] = ...
    MultiCellSchurInteriorSol(xSchur, Avvcell, InvAvvStorage, Avpcell, bvcell,...
    xCells, yCells, cellIndices)
    
    %% Description of the Input Arguments
    %  xSchur - reduced solution
    %  Avvcell - cell containing
    %  InvAvvStorage - cell containing LU factorizatoins of matrices in Avv
    %  Avpcell - couplings between Avv and App
    %  bvcell - source vector (Ax = b)
    %  xCell - number of metatoms in x direction
    %  yCells - number of metatoms in y direction
    %  SingleCellSize NxN nodes for each metatom
    %  cellindices
    
    tolSol = [xSchur];
    recSolIntArray = cell(xCells, yCells);
    %% Need a way to track all the correct cells to get the interior sol
    count = 1;
    for i = 1:xCells
        for j = 1:yCells
            Avv = Avvcell{count,count};
            L = InvAvvStorage{cellIndices(count), 1};
            U = InvAvvStorage{cellIndices(count), 2};
            bv = bvcell{count,1};
            Avp = Avpcell{count,1};
            yVec = bv- Avp*xSchur;
            recSolInt = (U\(L\yVec));
            tolSol = [tolSol; recSolInt];
            SingleCellSize = length(recSolInt)^.5;
            reshaped = reshape(recSolInt, SingleCellSize, SingleCellSize);
            recSolIntArray{i,j} = reshaped;
            count = count+1;
        end
    end
    
    
end