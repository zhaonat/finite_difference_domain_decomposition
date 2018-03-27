%% Matrix Operator Artificially Intelligent Row Swaps
% clear;
% close all;
% N(1) = 5; N(2) = 5;
% A2 = eye(N(1)^2,N(2)^2)
% A2(1,2) = 2;
% A2(2,1) = 2;

function [Q, SymA, symB,y, boundaryCells] = SingleCellBoundaryInterior(A2,b,N)
   
    indexorder = transpose(1:N(1)*N(2));
    
    %% Call IndexPermutation Function to achieve reordering
    [permutedIndices,interiorCellIndex, boundaryCellIndex, ...
        boundaryCells, interiorCells] = IndexPermutation(N);
    reversePerm = [boundaryCellIndex; interiorCellIndex];
    %% Execute a Symmetry Preserving Column Row Permutation Combination
    Q = speye(N(1)*N(2));
    for i = 1:N(1)*N(2)
       Q(i,i) = 0;
       indexshift = reversePerm(i);
       Q(i,indexshift) = 1;
    end

    %%test permutation matrix should permute index order to permuted indices
    y = Q*indexorder;

    %% Transform Equations
    SymA = Q*A2*transpose(Q);
    %issymmetric(SymA)
    symB = Q*b;
    
 
end