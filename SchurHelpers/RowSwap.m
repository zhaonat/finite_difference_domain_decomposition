%% Matrix Swapping Function
function A = RowSwap(A,i,j)

    %assert( i > 0 && i <= size(A,1) && j > 0 && j <= size(A,1) );
    A([i j],:) = A([j i],:);

end