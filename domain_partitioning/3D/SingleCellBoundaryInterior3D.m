%% Matrix Operator Artificially Intelligent Row Swaps
% clear;
% close all;
% N(1) = 5; N(2) = 5;
% A2 = eye(N(1)^2,N(2)^2)
% A2(1,2) = 2;
% A2(2,1) = 2;

function [y, Q, SymA, SymB, boundaryCells, interiorCells, hpart, vpart] = ...
    SingleCellBoundaryInterior3D(A,b,N)
    Nx = N(1); Ny = N(2); Nz = N(3);
    indexorder = transpose(1:N(1)*N(2)*N(3));
    
    %% Call IndexPermutation Function to achieve reordering
    [permutedIndices, boundaryCells, interiorCells] = IndexPermutation3D(N);
    %permutedIndicesSingle = IndexPermutation3D(N);
    hpart = length(boundaryCells);
    vpart = hpart;
    %% Execute a Symmetry Preserving Column Row Permutation Combination
    xind = zeros(3*Nx*Ny*Nz,1);
    yind = zeros(3*Nx*Ny*Nz,1);
    vals = ones(3*Nx*Ny*Nz,1);
    for i = 1:3*N(1)*N(2)*N(3)
       
       indexshift = permutedIndices(i);
       xind(i) = i;
       yind(i) = indexshift;
       %Q(i,indexshift) = 1;
    end
    Q = sparse(xind,yind,vals);
    
 
    SymA = Q*A*transpose(Q);
    SymB = Q*b;


%     Q2 = speye(N(1)*N(2)*N(3));
%     for i = 1:N(1)*N(2)*N(3)
%        Q2(i,i) = 0;
%        indexshift = permutedIndicesSingle(i);
%        Q2(i,indexshift) = 1;
%     end
%     Q3 = blkdiag(Q2,Q2,Q2);
%     %%test permutation matrix should permute index order to permuted indices
    
    indexorder = [indexorder;indexorder;indexorder];
    y = Q*indexorder;
    
    %% Transform Equations
   
%     
%     SymA2=Q3*A*transpose(Q3);
%     SymB2 = Q3*b;
%     
    

end