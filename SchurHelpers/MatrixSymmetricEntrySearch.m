%% Helper Function which finds and stores all indices and values for which a matrix is found to 
% be not symmetric
%A = [ 1 2 3; 4 5 6; 3 8 0];
function [mismatchDat grid coordinateMap] =  MatrixSymmetricEntrySearch(A)
    A = full(A);
    [xindices,yindices, values] = find(A) ;
    coordinateMap = [xindices, yindices, values];

    N = size(A);
    if(N(1)~=N(2))
       disp('matrix not square')
       mismatchDat = [];
       return
    end
    symmetry = true;
    mismatchx = []; mismatchy = []; 
    mismatchVal=[];

    for i = 1:N(1)
        for j = 1:N(2)
            if(abs(A(i,j)) ~= abs(A(j,i)))
               mismatchx = [mismatchx i];
               mismatchy = [mismatchy j];
               mismatchVal = [mismatchVal; A(i,j) A(j,i)];
            end
        end
    end

    mismatchDat = [transpose(mismatchx) transpose(mismatchy) mismatchVal];

    %% Return a grid which visualizes precise locations of mismatch
    grid = zeros(N(1), N(2));
    for i = 1:length(mismatchx)
        x = mismatchx(i);
        y = mismatchy(i);
        grid(x,y) = 1;
    end

end