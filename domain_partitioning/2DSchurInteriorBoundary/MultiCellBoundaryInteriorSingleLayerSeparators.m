%% MulticellRegime with one partition for the boundary, one for the interio nodes

%for now, we only consider square cell subdivisions for simplicity
function [SymA, symB, Q, permutedIndices, boundaryCells, interiorCells, hpart, vpart] = ...
    MultiCellBoundaryInteriorSingleLayerSeparators(A, b, xCells, yCells,SingleCellSize, N)
M = N(1)*N(2);
indexorder = transpose(1:N(1)*N(2));
NaturalOrdering = CoordinateIndexing(N(1), N(2));

%figure;
%imagesc(A2)


Nx = N(1); Ny = N(2);
boundaryIndices = [1];
for i = 1:xCells
    boundaryIndices = [boundaryIndices, boundaryIndices(i)+SingleCellSize+1];
end

boundaryCells = [];
interiorCells = [];
boundaryCellIndex = []; interiorCellIndex = [];

for j = 1:xCells
    for k = 1:yCells
        xend = 1+SingleCellSize*j+j
        yend = 1+SingleCellSize*k+k
        xbeg = 1+SingleCellSize*(j-1)+(j-1)
        ybeg = 1+SingleCellSize*(k-1)+(k-1)
        if(xbeg > 1)
           xbeg = xbeg+1; 
        end
        if(ybeg > 1)
            ybeg = ybeg+1
        end
        for i = 1:(Nx*Ny)
           startIndex = NaturalOrdering(i,:); 
           x = startIndex(1); y = startIndex(2);
           if(x>=xbeg && x<= xend && y >= ybeg && y <= yend)


               if(any(x == boundaryIndices)||any(y==boundaryIndices))
                   boundaryCells = [boundaryCells; startIndex];
                   boundaryCellIndex = [boundaryCellIndex; i];

               else

                   interiorCells = [interiorCells; startIndex];
                   interiorCellIndex = [interiorCellIndex; i];
               end
           else
               continue
           end
        end
    
    end
end

vpart = length(boundaryCellIndex);
hpart = vpart;
permutedIndices = [boundaryCellIndex; interiorCellIndex];

%%Permute Indices one more time so that all the interior cells are grouped
%%correctly

%% Execute a Symmetry Preserving Column Row Permutation Combination
 xind = zeros(Nx*Ny,1);
yind = zeros(Nx*Ny,1);
vals = ones(Nx*Ny,1);

%% ALWAYS CREATE THE PERMUTATION SPARSE MATRIX LIKE THIS!! 
for i = 1:N(1)*N(2)
   indexshift = permutedIndices(i);
   xind(i) = i;
   yind(i) = indexshift;
   %Q(i,indexshift) = 1;
end
Q = sparse(xind,yind,vals);

%% Transform Equations
SymA = Q*A*transpose(Q);
%issymmetric(SymA)
symB = Q*b;

end