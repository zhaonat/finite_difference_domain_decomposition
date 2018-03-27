%% 3D Schur Complement Reordering Function

function [permutedIndices, boundaryCells, interiorCells, hpart, vpart, l1, l2] =...
    IndexPermutation3DTwoCell(N,div)
    %N = [10 10 10];
    NaturalOrdering = CoordinateIndexing3D(N(1), N(2), N(3));
    Nx = N(1); Ny = N(2); Nz = N(3);
    boundaryCellIndex = [];
    interiorCellIndex = []; interiorCellIndex2 = [];
    edgeCellIndex = [];
    
    boundaryCells = [];
    interiorCells = []; interiorCells2 = [];
    inter = 1; bounder = 1; bounderDiv = 1;
    M= 3*N(1)*N(2)*N(3);
    Vol = N(1)*N(2)*N(3);
    for i = 1:3*(N(1)*N(2)*N(3))
       
       ind = mod(i,M/3);
       if(ind == 0)
           ind = M/3;
       end
       %% this differentiates which E-field component the index belongs too
       startIndex = NaturalOrdering(ind,:); 
       if(i <= Vol)
           startIndex(1) =startIndex(1)+ 1/2;
       elseif(i> Vol && i<= 2*Vol)
           startIndex(2)= startIndex(2) + 1/2;
       else
           startIndex(3) = startIndex(3) + 1/2;
       end
       
       endpoints = max(NaturalOrdering);        
       xend = endpoints(1); yend = endpoints(2); zend = endpoints(3);
       x = startIndex(1); y = startIndex(2); z = startIndex(3);

       if( x == 1.5 || x == xend+0.5 || y == 1.5 || y == yend+0.5 || ...
               z == 1.5 || z == zend+0.5 ||x == 1 || y == 1 || z == 1 || x == xend ...
               || z == zend || y == yend)

           boundaryCells=[ boundaryCells; startIndex];
           boundaryCellIndex = [boundaryCellIndex; i];
           bounder = bounder+1;
       
       %% create a partition along the x-axis;
       elseif( x == Nx/div || x == Nx/div + 0.5 || x == Nx/div+1 )%|| x == Nx/div -0.5)
           boundaryCells=[ boundaryCells; startIndex];
           boundaryCellIndex = [boundaryCellIndex; i];
           bounderDiv = bounderDiv+1;    
       else
           if(x <= Nx/div-0.5)
               interiorCells = [interiorCells; startIndex];
               interiorCellIndex = [interiorCellIndex; i];
           else
               interiorCells2 = [interiorCells2; startIndex];
               interiorCellIndex2 = [interiorCellIndex2; i]; 
           end
           inter = inter+1;
       end
    end
    l1 = length(interiorCells);
    l2 = length(interiorCells2);
    interiorCells = [interiorCells; interiorCells2];
    permutedIndices = [boundaryCellIndex; interiorCellIndex; interiorCellIndex2];
    hpart = length(boundaryCellIndex);
    vpart = hpart;
 end