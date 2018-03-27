%% Triple Field Permutation

function [permutedIndices, boundaryCells, interiorCells] = IndexPermutation3D_3Field(N)
    %N = [10 10 10];
    NaturalOrdering = CoordinateIndexing3D(N(1), N(2), N(3));
    M = 3*N(1)*N(2)*N(3);


    boundaryCellIndex = [];
    interiorCellIndex = [];
    edgeCellIndex = [];
    
    boundaryCells = [];
    interiorCells = [];
    inter = 1; bounder = 1;
    Vol = N(1)*N(2)*N(3);
    for i = 1:3*(N(1)*N(2)*N(3)) %we will move the boundaries for all three fields to the end
       ind = mod(i,M/3);
       if(ind == 0)
           ind = M/3;
       end
       startIndex = NaturalOrdering(ind,:); 
       x = startIndex(1); y = startIndex(2); z = startIndex(3);
       if(i <= Vol)
           startIndex(1) =startIndex(1)+ 1/2;
       elseif(i> Vol && i<= 2*Vol)
           startIndex(2)= startIndex(2) + 1/2;
       else
           startIndex(3) = startIndex(3) + 1/2;
       end

       if(i <= Vol && (z == 1 || x==N(1) || y == 1 ))
           
           boundaryCells=[ boundaryCells; startIndex];
           boundaryCellIndex = [boundaryCellIndex; i];
           bounder = bounder+1;
%        elseif(mod(x, N(1)) == N(1)-1 || mod(x,N(1)) == 2 || mod(y, N(2)) == N(2)-1 || mod(y,N(2)) == 2 || ...
%                mod(z, N(3)) == N(3)-1 || mod(z,N(3)) == 2 )
%            edgeCellIndex = [edgeCellIndex; i];
       elseif( i> Vol && i<= 2*Vol &&( y ==N(2) || x == 1 || z== 1 ))
          
           boundaryCells=[ boundaryCells; startIndex];
           boundaryCellIndex = [boundaryCellIndex; i];
           bounder = bounder+1;
       elseif( i> 2*Vol && i<= 3*Vol &&( z ==N(3) || x == 1 || y== 1 ))
         
           boundaryCells=[ boundaryCells; startIndex];
           boundaryCellIndex = [boundaryCellIndex; i];
           bounder = bounder+1;
       else
           interiorCells = [interiorCells; startIndex];
           interiorCellIndex = [interiorCellIndex; i];
           inter = inter+1;
       end
    end
    permutedIndices = [boundaryCellIndex; edgeCellIndex; interiorCellIndex];
    
 end