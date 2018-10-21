%% Matrix Operator Artificially Intelligent Row Swaps
% clear;
% close all;
% N(1) = 5; N(2) = 5;
% A2 = eye(N(1)^2,N(2)^2)
% A2(1,2) = 2;
% A2(2,1) = 2;

%for now, we only consider square cell subdivisions for simplicity
function [SymA, symB, Q, permutedIndices, boundaryCells, InteriorCells] = ...
    MultiCellBoundaryInterior(A,b, divx, divy,N)
M = N(1)*N(2);
indexorder = transpose(1:N(1)*N(2));
NaturalOrdering = CoordinateIndexing(N(1), N(2));

%figure;
%imagesc(A2)

modValuex = N(1)/divx;
modValuey = N(2)/divy;

boundaryCells = [];
interiorCells = [];
boundaryCellIndex = []; interiorCellIndex = [];
for i = 1:ceil(N(1)*N(2))
   startIndex = NaturalOrdering(i,:); 
   if(i == 1 || i <= N(1) || mod(i, N(1)) == 0 || mod(i,N(1)) == 1 || (i>= M-N(1)+1 && i<=M))
       boundaryCells = [boundaryCells; startIndex];
       boundaryCellIndex = [boundaryCellIndex; i];
       
   elseif(mod(i, modValuex) == 0 || mod(i, modValuex) == 1)
       boundaryCells = [boundaryCells; startIndex];
       boundaryCellIndex = [boundaryCellIndex; i];
   elseif(mod(i, modValuey) == 0 || mod(i, modValuey) == 1)
       boundaryCells = [boundaryCells; startIndex];
       boundaryCellIndex = [boundaryCellIndex; i];
   else
       
       interiorCells = [interiorCells; startIndex];
       interiorCellIndex = [interiorCellIndex; i];
   end
end

permutedIndices = [boundaryCellIndex; interiorCellIndex];

%% Execute a Symmetry Preserving Column Row Permutation Combination
Q = speye(N(1)*N(2));
for i = 1:N(1)*N(2)
   Q(i,i) = 0;
   indexshift = permutedIndices(i);
   Q(i,indexshift) = 1;
end

%%test permutation matrix should permute index order to permuted indices
y = Q*indexorder;

%% Transform Equations
SymA = Q*A*transpose(Q);
%issymmetric(SymA)
symB = Q*b;

end