function [Nindexing, Rindexing] = CoordinateIndexing(N1, N2)

    Nindexing = zeros(N1*N2,2);
    Rindexing = zeros(N1, N2);
    counter = 1;
    for i = 1:N1;
       x = i;
       for j = 1:N2;
           y = j;
           point = [x,y];
           Nindexing(counter,:) =point ;
           Rindexing(x,y) = counter;
           counter = counter+1;
       end
    end

end