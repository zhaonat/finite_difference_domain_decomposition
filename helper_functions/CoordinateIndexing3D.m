function [IndexToCoord, CoordToIndex] = CoordinateIndexing3D(N1, N2, N3)

    IndexToCoord = zeros(N1*N2*N3,3);
    CoordToIndex = zeros(N1, N2, N3);
    counter = 1;
    for i = 1:N1;
       x = i;
       for j = 1:N2;
           y = j;
           for k = 1:N3
                z = k;
                point = [x, y, z];
                IndexToCoord(counter,:) = point;
                CoordToIndex(x,y,z) = counter;
                counter = counter+1;
           end
       end
    end

end