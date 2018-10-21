%% Single layer Schur complement with PML setup
close all
clear

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.55*L0;  % wavelength in L0
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 40;
epsilon = 10;
Sx = SingleCellSize; Sy = SingleCellSize;

k = 2; xCells =k; yCells = k;

%% ================ Domain Size Parameters ==========================
N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
N0 = N;
Npml = [10 10];  
xrange = k*L0*[-1 1];  % x boundaries in L0
yrange = k*L0*[-1 1];  % y boundaries in L0

%% Generate Parameters for Domain with PML
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  

%% ==============Setup Dielectric ===================================
[eps_air, cellIndices] = ... 
    multiRandomCellDielectricSingleLayerSep(k, k, SingleCellSize,...
    SingleCellSize, Npml,epsilon); %% ADD coe to account for PML

%% setup point source 
Mz = zeros(N);
ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

Mz1 = zeros(N0);
Mz1(ceil(N0/2)) = 1;
%% ============ ==== FDFD Matrix Setup ==== ===========================%
%with pml
[A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices_unscaled(wvlen, xrange, ...
    yrange, eps_air, Mz, Npml);


%% ===========TEST OF THE REORDERING ================================%

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
cx = 1; cy = 1;
for i = 1:(Nx*Ny) %THIS IS HIGHLY INEFFICIENT
   startIndex = NaturalOrdering(i,:); 
   x = startIndex(1); y = startIndex(2);
   key=strcat(num2str(x),',', num2str(y));

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
pmlxpart = M - length(pmlCellIndex);
pmlypart = pmlxpart;

permutedIndices = [boundaryCellIndex; interiorCellIndex; pmlCellIndex];

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



%% ============= Test Solve ==============================
x2 = A\b;
Hz = reshape(x2, length(x2)^.5, length(x2)^.5);
imagesc(abs(Hz))
hold on
scatter(pmlCell(:,1), pmlCell(:,2))
scatter(boundaryCells(:,1), boundaryCells(:,2))
scatter(interiorCells(1:1600,1), interiorCells(1:1600,2), '.')

%% Visualize the ordering of the pml:
figure()
scatter(boundaryCells(:,1), boundaryCells(:,2))
hold on
scatter(pmlCell(1:800,1), pmlCell(1:800,2), '.')


