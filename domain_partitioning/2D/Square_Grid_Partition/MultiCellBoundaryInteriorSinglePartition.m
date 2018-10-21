%% MulticellRegime with one partition for the boundary, one for the interio nodes
%% There is no PML (only pbc or dirichlet or normal boundary values)
%for now, we only consider square cell subdivisions for simplicity
function [Q, SymA, symB, permutedIndices, boundaryCells, interiorCells, hpart, vpart] = ...
    MultiCellBoundaryInteriorSinglePartition(A, b, xCells, yCells,SingleCellSize, N)

    Npml = [0,0];
        %% first, let's get back the dimensions of the small domain
    Noriginal = N-2*Npml; %factor of 2 accounts for the fact that you need a pml on each side
    pmloffset = 2*Npml(1);
    pmlwidth = Npml(1);
    M = N(1)*N(2);
    indexorder = transpose(1:N(1)*N(2));
    NaturalOrdering = CoordinateIndexing(N(1), N(2));

    Nx = N(1); Ny = N(2);
    boundaryIndices = [1];
    for i = 1:xCells
        boundaryIndices = [boundaryIndices, boundaryIndices(i)+...
            SingleCellSize+1];
    end
    boundaryIndices = pmloffset/2 + boundaryIndices;

    Cx = length(boundaryIndices)-1;
    totCells = Cx^2;
    %% these cells store the actual x and y coordinates
    boundaryCells = [];
    interiorCellIndexDict = cell(totCells,1); 
    interiorCellDict = cell(totCells, 1);

    pmlCell = [];
    boundaryCellIndex = []; interiorCellIndex = [];
    pmlCellIndex = [];

    %% we need a map of all the boundaryindices in x and y to cell index
    boundaryMap = []; counter = 1;
    for i = boundaryIndices(1:end-1)
        for j = boundaryIndices(1:end-1)
            boundaryMap = [boundaryMap; i,j, counter];
            counter = counter+1
        end
    end
    nodeMap = zeros(N);
    d = size(boundaryMap);
    for i = 1:d(1)
        nodeMap(boundaryMap(i,1)+1:boundaryMap(i,1)+SingleCellSize,...
            boundaryMap(i,2)+1:boundaryMap(i,2)+SingleCellSize) = ...
        boundaryMap(i,3);
    end

    %% cycle through all points on the grid in integer indexing
    cellIndex = 1;
    for i = 1:(Nx*Ny) %THIS IS HIGHLY INEFFICIENT, ideally want something that's vectorized
       startIndex = NaturalOrdering(i,:); 
       x = startIndex(1); y = startIndex(2);

       %% condition for being inside the pml buffer
       if(x > pmlwidth+Noriginal(1) || y > pmlwidth+Noriginal(2) ...
           || x <= Npml(1) || y <= Npml(2))
           pmlCell = [pmlCell; startIndex];
           pmlCellIndex = [pmlCellIndex; i];

       %% condition for being inside the structure
       else
           if(any(x == boundaryIndices)||any(y==boundaryIndices))
               boundaryCells = [boundaryCells; startIndex];
               boundaryCellIndex = [boundaryCellIndex; i];

           else
               cellIndex = nodeMap(x,y);
               interiorCellDict{cellIndex} = [interiorCellDict{cellIndex}; startIndex];
               interiorCellIndexDict{cellIndex}=...
                   [interiorCellIndexDict{cellIndex}; i];
           end

       end

    end %% end of natural ordering indexing for loop

    interiorCells = cell2mat(interiorCellDict);
    interiorCellIndex = cell2mat(interiorCellIndexDict);


    vpart = length(boundaryCellIndex);
    hpart = vpart;
    permutedIndices = [boundaryCellIndex; interiorCellIndex; pmlCellIndex];
    pmlxpart = M - length(pmlCellIndex);
    pmlypart = pmlxpart
    %%Permute Indices one more time so that all the interior cells are grouped
    %%correctly

    %% Execute a Symmetry Preserving Column Row Permutation Combination
    xind = zeros(Nx*Ny,1);
    yind = zeros(Nx*Ny,1);
    vals = ones(Nx*Ny,1);

    %% ALWAYS CREATE THE PERMUTATION SPARSE MATRIX LIKE THIS!! 
    for i = 1:N(1)*N(2)
       indexshift = permutedIndices(i);
       xind(i) = i;
       yind(i) = indexshift;
       %Q(i,indexshift) = 1;
    end
    Q = sparse(xind,yind,vals);

    %% Transform Equations
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    symB = Q*b;

end