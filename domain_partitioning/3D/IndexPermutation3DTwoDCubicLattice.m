
function [permutedIndices, boundaryCells, interiorCells, hpart, vpart] = IndexPermutation3DTwoDCubicLattice(N,divx, divy)
    %N = [10 10 10];
    NaturalOrdering = CoordinateIndexing3D(N(1), N(2), N(3));
    Nx = N(1); Ny = N(2); Nz = N(3);
    boundaryCellIndex = [];
    interiorCellIndex = [];
    edgeCellIndex = [];
    
    boundaryCells = [];
    interiorCells = [];
    inter = 1; bounder = 1;
    M= 3*N(1)*N(2)*N(3);
    Vol = N(1)*N(2)*N(3);
    for i = 1:3*(N(1)*N(2)*N(3))
       
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
       elseif( x == Nx/divx || x == Nx/divx -0.5 || x == Nx/divx+0.5 || x == Nx/divx+1)
           boundaryCells=[ boundaryCells; startIndex];
           boundaryCellIndex = [boundaryCellIndex; i];
           bounder = bounder+1;    
       elseif( y == Ny/divy || y == Ny/divy -0.5 || y == Ny/divy+0.5 || y == Ny/divy+1)
           boundaryCells=[ boundaryCells; startIndex];
           boundaryCellIndex = [boundaryCellIndex; i];
           bounder = bounder+1;    
       else
           interiorCells = [interiorCells; startIndex];
           interiorCellIndex = [interiorCellIndex; i];
           inter = inter+1;
       end
    end
    permutedIndices = [boundaryCellIndex; interiorCellIndex];
    hpart = length(boundaryCellIndex);
    vpart = hpart;
 end