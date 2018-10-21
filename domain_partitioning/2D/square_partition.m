
function [Q,mask] = square_partition(N, coord1, coord2, exist_mask)
    %N = [Nx,Ny]
    % coord1 = (x1, y1)
    % coord2 = (x2, y2); x2 > x1, y2 > y1
    % also return a mask, showing where the subdomains are, 
    % this function is not additive...meaning we can't do this if there 
    % are multiple disconnected squares on our grid
    
    %cell containing the indices of all the 
    
    x1 = coord1(1); y1 = coord1(2);
    x2 = coord2(1); y2 = coord2(2);
    
    [xx,yy] = meshgrid(1:N(1), 1:N(2));
    M = prod(N);
    %natural_order_grid = reshape(1:M, N(1), N(2));
    
    [subdomain_mask_xx, subdomain_mask_yy] = meshgrid(1:N(1), 1:N(2));
    
    subdomain_mask_xx(xx<x2 & xx >x1 &yy > y1 & yy < y2) = 0;
    subdomain_mask_yy(yy > y1 & yy < y2 & xx<x2 & xx >x1) = 0;
    subdomain_mask_yy( subdomain_mask_yy >0) = 1;
    mask = subdomain_mask_yy+exist_mask;
    mask(mask <2) = 0;
    mask(mask == 2) = 1;

    interior = find(subdomain_mask_xx == 0| exist_mask == 0);
    exterior = find(subdomain_mask_xx>0 & exist_mask >0);
%    permutation_cell = cell(1);
%     permutation_cell{1} = interior;
%     permutation_cell{2} = exterior;
    permutation = [interior; exterior];
    
    %% construct permutation matrix;
    Q = sparse(1:M,permutation.',ones(M,1));
    
    %% recall the correct permutation operation is:
%     SymA = Q*A*transpose(Q);
%     symB = Q*b;
    
end