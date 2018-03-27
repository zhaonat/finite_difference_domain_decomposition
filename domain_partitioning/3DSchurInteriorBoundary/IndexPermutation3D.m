%% 3D Schur Complement Reordering Function

function [permutedIndices, boundaryCells, interiorCells] = IndexPermutation3D(N)
    %N = [10 10 10];
    NaturalOrdering = CoordinateIndexing3D(N(1), N(2), N(3));
    Nx = N(1); Ny = N(2); Nz = N(3);
    boundaryCellIndex = zeros(3*6*Nx*Ny,1);
    interiorCellIndex = zeros(3*Nx*Ny*Nz,1);
    endpoints = max(NaturalOrdering);        
    xend = endpoints(1); yend = endpoints(2); zend = endpoints(3);
    boundaryCells = zeros(3*6*Nx*Ny,3);
    interiorCells = zeros(3*Nx*Ny*Nz,3);
    inter = 1; bounder = 1;
    M= 3*N(1)*N(2)*N(3);
    Vol = N(1)*N(2)*N(3);
    for i = 1:3*(N(1)*N(2)*N(3))
       %% this order should successfully cycle through Ex, Ey, then Ez
       ind = mod(i,M/3);
       if(ind == 0)
           ind = M/3;
       end
       startIndex = NaturalOrdering(ind,:); 
       if(i <= Vol)
           startIndex(1) =startIndex(1)+ 1/2;
       elseif(i> Vol && i<= 2*Vol)
           startIndex(2)= startIndex(2) + 1/2;
       else
           startIndex(3) = startIndex(3) + 1/2;
       end
       
       x = startIndex(1); y = startIndex(2); z = startIndex(3);
         
       if( x == 1.5 || x == xend+0.5 || y == 1.5 || y == yend+0.5 || ...
               z == 1.5 || z == zend+0.5 ||x == 1 || y == 1 || z == 1 || x == xend ...
               || z == zend || y == yend)

           boundaryCells(bounder,1:3)= startIndex;
           boundaryCellIndex(bounder,1) = i;
           bounder = bounder+1;

       else
           interiorCells(inter,1:3) = startIndex;
           interiorCellIndex(inter,1) = i;
           inter = inter+1;
       end
    end
    boundaryCellIndex = boundaryCellIndex(boundaryCellIndex~=0);
    interiorCellIndex = interiorCellIndex(interiorCellIndex~=0);
    boundaryCells( all(~boundaryCells,2), : ) = [];
    interiorCells( all(~interiorCells,2), : ) = [];
   permutedIndices = [boundaryCellIndex; interiorCellIndex];
    
 end