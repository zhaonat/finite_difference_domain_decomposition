%% 3D Schur Complement Reordering Function

function [permutedIndices, boundaryCells, interiorCells,interiorCellStorage, hpart, vpart ] =...
    IndexPermutation3DOneDChain(N,divx)
    %N = [10 10 10];
    NaturalOrdering = CoordinateIndexing3D(N(1), N(2), N(3));
    Nx = N(1); Ny = N(2); Nz = N(3);
    boundaryCellIndex = [];
    interiorCellIndex = [];
    interiorCellStorage = cell(divx,1);
    
    boundaryCells = [];
    interiorCells = [];
    inter = 1; bounder = 1; bounderDiv = 1;
    M= 3*N(1)*N(2)*N(3);
    Vol = N(1)*N(2)*N(3);
    xLen = Nx/divx; xCells = divx; 
    yLen = Ny; zLen = Nz;
    for k = 0:xCells-1

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
    
       if(x <=(k+1)*xLen+0.5 && x > k*xLen+0.5)
        
           
            if(((mod(x, xLen) == 0 ||  mod(x-0.5,xLen) == 0 ...
                    || mod(x-1,xLen) ==0 || mod(x-1,xLen)== 0)  ||...
           ( x == 1.5 || x == xend+0.5 || y == 1.5 || y == yend+0.5 || ...
                   z == 1.5 || z == zend+0.5 ||x == 1 || y == 1 || z == 1 || x == xend ...
                   || z == zend || y == yend|| x == 2 || x == 2.5))) %get rid of this and you get just x direction partitions

               boundaryCells=[ boundaryCells; startIndex];
               boundaryCellIndex = [boundaryCellIndex; i];
               bounder = bounder+1;

            elseif(y > 1.5 && y < yend  && ...
                   z> 1.5 && z < zend )
               interiorCells = [interiorCells; startIndex];
               interiorCellIndex = [interiorCellIndex; i];
               interiorCellStorage{k+1,1} = [interiorCellStorage{k+1,1},i];
               inter = inter+1;
               
            end
               
        else
            continue
        end
           
    
       end
       
       
     end
    
    permutedIndices = [boundaryCellIndex; interiorCellIndex];
    hpart = length(boundaryCellIndex);

    vpart = hpart;
 end