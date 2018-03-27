
function [SymA,SymB, Q, permutedIndices, interiorCoord, boundaryCoord, hpart, vpart] = ...
MultiCellReorder(A, b,xCells, yCells, N)

    Nx = N(1); 
    Ny = N(2);

    xDiv = Nx/xCells;
    yDiv = Ny/yCells;
    totCells = xCells*yCells;
    %NaturalOrdering = CoordinateIndexing(Nx, Ny); %% expensive? not really.
    boundaryInd = zeros(totCells*(2*xDiv+2*yDiv-4),1); 
    interiorInd = zeros(Nx*Ny -totCells*(2*xDiv+2*yDiv-4),1);
    boundaryCoord = zeros(totCells*(2*xDiv+2*yDiv-4),2); 
    interiorCoord = zeros(Nx*Ny -totCells*(2*xDiv+2*yDiv-4),2);
    bc = 1; ic = 1;
    Q = speye(N(1)*N(2));

    for i = 0:xCells-1
       for j = 0:yCells-1
            xBeg = i*xDiv+1;
            yBeg = j*yDiv+1;
            xEnd = (i+1)*xDiv;
            yEnd = (j+1)*yDiv;
            counter  = 1;
            for k = 1:Nx*Ny
               %Q(k,k) = 0;
               [x,y] = NatIndexToCoord(k,Nx, Ny);
%                x = NaturalOrdering(k,1);
%                y = NaturalOrdering(k,2);
               if(x > xBeg && x< xEnd && y> yBeg && y < yEnd)
                  interiorInd(ic) = k; 
                  interiorCoord(ic,:) = [x,y];
                  %Q(counter,k) = 1;
                  ic = ic+1;
               elseif((x == xBeg && y >= yBeg && y <=yEnd) || (x == xEnd && y >= yBeg&& y <=yEnd) || ...
                       (y == yBeg && x >= xBeg && x <= xEnd) || (y == yEnd && x >= xBeg && x <= xEnd))
                   boundaryInd(bc) = k;
                   %Q(counter,k) = 1;
                   boundaryCoord(bc,:) = [x,y];
                   bc = bc+1;
                   
               end
               counter = counter+1;
            end
       end
    end
    permutedIndices = [boundaryInd; interiorInd];
    
    %% Get the new matrix Q
    xind = zeros(Nx*Ny,1);
    yind = zeros(Nx*Ny,1);
    vals = ones(Nx*Ny,1);
    for i = 1:N(1)*N(2)
       
       indexshift = permutedIndices(i);
       xind(i) = i;
       yind(i) = indexshift;
       %Q(i,indexshift) = 1;
    end
    Q = sparse(xind,yind,vals);
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    SymB = Q*b;
    
    hpart = length(boundaryInd);
    vpart = hpart;

    
end