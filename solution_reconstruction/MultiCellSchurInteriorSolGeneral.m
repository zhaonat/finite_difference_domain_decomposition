%% THIS FUNCTION CAN BE USED FOR BOTH SINGLE AND DUAL LAYER SEPARATORS!
%% the interior solution really just needs InvAvvStorage

function tolSol = ...
    MultiCellSchurInteriorSolGeneral(xSchur, Avvcell, ...
    Avpcell, bvcell)
    
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
    
    tolSol = [xSchur];
    %% Need a way to track all the correct cells to get the interior sol
    for count = 1:length(bvcell)
        Avv = Avvcell{count,count};
        bv = bvcell{count};
        Avp = Avpcell{count};
        yVec = bv - Avp*xSchur;
%         L = InvAvvStorage{count, 1};
%         U = InvAvvStorage{count, 2};
        recSolInt = Avv\yVec;
        tolSol = [tolSol; recSolInt];
        %this explicitly assuemes that all domains we Schur
        %complemented are square in shape...

    end
    
    
end