N = [10 10 10];
    NaturalOrdering = CoordinateIndexing3D(N(1), N(2), N(3));
    M = 3*N(1)*N(2)*N(3);


    boundaryCellIndex = [];
    interiorCellIndex = [];
    boundaryCells = [];
    interiorCells = [];
    inter = 1; bounder = 1;
    for i = 1:3*(N(1)*N(2)*N(3))
       startIndex = NaturalOrdering(mod(i, M/3)+1,:); 
       x = startIndex(1); y = startIndex(2); z = startIndex(3);

       if( mod(x, N(1)) == 0 || mod(x,N(1)) == 1 || mod(y, N(2)) == 0 || mod(y,N(2)) == 1 || ...
               mod(z, N(3)) == 0 || mod(z,N(3)) == 1  )

           boundaryCells=[ boundaryCells; startIndex];
           boundaryCellIndex = [boundaryCellIndex; i];
           bounder = bounder+1;
       else
           interiorCells = [interiorCells; startIndex];
           interiorCellIndex = [interiorCellIndex; i];
           inter = inter+1;
       end
    end
    permutedIndices = [boundaryCellIndex; interiorCellIndex];