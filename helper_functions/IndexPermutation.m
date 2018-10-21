
function [permutedIndices, interiorCellIndex, boundaryCellIndex, boundaryCells, interiorCells] = IndexPermutation(N)
    NaturalOrdering = CoordinateIndexing(N(1), N(2));
    M = N(1)*N(2);
    %figure;
    %imagesc(A2)
    boundaryCount = 2*N(1)+2*(N(2)-2);
    interiorCount = M-boundaryCount;
    boundaryCellIndex = zeros(boundaryCount,1);
    interiorCellIndex = zeros(interiorCount,1);
    boundaryCells = zeros(boundaryCount,2);
    interiorCells = zeros(interiorCount,2);
    inter = 1; bounder = 1;
    for i = 1:ceil(N(1)*N(2))
       startIndex = NaturalOrdering(i,:); 
       if(i == 1 || i <= N(1) || mod(i, N(1)) == 0 || mod(i,N(1)) == 1 || (i>= M-N(1)+1 && i<=M))
           boundaryCells(bounder,:) = startIndex;
           boundaryCellIndex(bounder) = i;
           bounder = bounder+1;
       else
           interiorCells(inter,:) = startIndex;
           interiorCellIndex(inter) = i;
           inter = inter+1;
       end
    end
    permutedIndices = [interiorCellIndex; boundaryCellIndex];

end