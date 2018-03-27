function [SymA, SymB, Q, permutedIndices, hpart, vpart, map] =...
    RectangularReorder(A,b, N, bounds, order) 
  
    %% order is a parameter which dictates direction of complement
    % order = 0 (interior to exterior)
    % order = 1 (exterior to interior)
    
    % bounds should be organized as
    %[xmin, ymin, xmax, ymax];

    %% get the coordinate indexing
    [~, OrderMap] = CoordinateIndexing(N(1), N(2));
    Nx = N(1); Ny = N(2);
    
    %% visualize Schur complement
    map = zeros(Nx, Nx);
    map(bounds(1):bounds(3), bounds(2):bounds(4)) = 1;

    %% Create Storage Data Structures
    interior = OrderMap(bounds(1):bounds(3), bounds(2):bounds(4));
    exterior1 = OrderMap(1:bounds(1)-1,:);
    exterior2 = OrderMap(bounds(3)+1:end, :);
    exterior3 = OrderMap(bounds(1):bounds(3), 1:bounds(2)-1);
    exterior4 = OrderMap(bounds(1):bounds(3), bounds(4)+1:end);
    map(1:bounds(1)-1,:) = 2;
    map(bounds(3)+1:end, :) = 3;
    map(bounds(1):bounds(3), 1:bounds(2)-1) = 4;
    map(bounds(1):bounds(3), bounds(4)+1:end) = 5;
    

    %% reshape exteriors
    exterior1 = reshape(exterior1, numel(exterior1),1);
    exterior2 = reshape(exterior2, numel(exterior2),1);
    exterior3 = reshape(exterior3, numel(exterior3),1);
    exterior4 = reshape(exterior4, numel(exterior4),1);
    
    interiorCellIndex = reshape(interior, numel(interior),1);
    exteriorCellIndex = [exterior1;exterior2; exterior3; exterior4];
    
    if(order == 1)
        permutedIndices = [interiorCellIndex; exteriorCellIndex];
        %% Partition Sizes
        hpart = length(interiorCellIndex);
        vpart = hpart;
    else
        permutedIndices = [exteriorCellIndex; interiorCellIndex];
        hpart = length(exteriorCellIndex);
        vpart = hpart;
    end

    
    %% Execute a Symmetry Preserving Column Row Permutation Combination
    xind = zeros(Nx*Ny,1);
    yind = zeros(Nx*Ny,1);
    vals = ones(Nx*Ny,1);

    %% ALWAYS CREATE THE PERMUTATION SPARSE MATRIX LIKE THIS!! 
    for i = 1:N(1)*N(2)
       indexshift = permutedIndices(i);
       xind(i) = i;
       yind(i) = indexshift;
    end
    
    Q = sparse(xind,yind,vals);

    %% Transform Equations
    size(Q)
    size(A)
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    SymB = Q*b; 
    
    
end