%% function to reorder indices in the multicell regime so that the ordering
%% lists indices in order of membership of cell

N = [100 100];
xCells = 2; yCells = 2;

    Nx = N(1); 
    Ny = N(2);

    xDiv = Nx/xCells;
    yDiv = Ny/yCells;

    NaturalOrdering = CoordinateIndexing(N(1), N(2));
    permutedIndices = [];
    boundaryInd = []; interiorInd = []; boundaryCoords = [];
    
    for i = 0:xCells-1
       for j = 0:yCells-1
            xBeg = i*xDiv+1
            yBeg = j*yDiv+1
            xEnd = (i+1)*xDiv
            yEnd = (j+1)*yDiv

            for k = 1:length(NaturalOrdering)
               x = NaturalOrdering(k,1);
               y = NaturalOrdering(k,2);
               if(x > xBeg && x< xEnd && y> yBeg && y < yEnd)
                  interiorInd = [interiorInd k]; 
               elseif((x == xBeg && y >= yBeg && y <=yEnd) || (x == xEnd && y >= yBeg&& y <=yEnd) || ...
                       (y == yBeg && x >= xBeg && x <= xEnd) || (y == yEnd && x >= xBeg && x <= xEnd))
                   boundaryInd = [boundaryInd k];
                   boundaryCoords = [boundaryCoords; x, y];
               end
            end
       end
    end
    permutedIndices = [boundaryInd, interiorInd];
    
    %% Get the new matrix Q
    Q = speye(N(1)*N(2));
    for i = 1:N(1)*N(2)
       Q(i,i) = 0;
       indexshift = permutedIndices(i);
       Q(i,indexshift) = 1;
    end
    
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    SymB = Q*b;
    
    hpart = length(boundaryInd);
    vpart = hpart;
    
    Grid = zeros(Nx, Ny);
    for i = 1:length(boundaryCoords)
        xi =boundaryCoords(i,1); 
        yi = boundaryCoords(i,2);
       Grid(xi, yi) = Grid(xi,yi)+ 1; 
    end
    
    imagesc(Grid)
