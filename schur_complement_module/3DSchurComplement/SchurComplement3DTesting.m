close all
clear;

%% Parameter Set up
%% Parameter Set-up
L0 = 1e-6;  % length unit: microns
wvlen = 4.0;  % wavelength in L0


for cellDim = 2:2
    cellx = cellDim; celly = cellDim; cellz = cellDim;

    Npml = [0 0 0];  % [Nx_pml Ny_pml]
    xrange = [-1 1];  % x boundaries in L0
    yrange = [-1 1];  % y boundaries in L0
    zrange = [-1 1];
    %[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  
    % domain is expanded to include PML

    %% Set up the permittivity.
    diel = 3;
    inclusionSize = 1;
    SingleCellDims = [5,5,5];
    [eps_r, interiorCoords, borderCoords] = cubeDielectricGrid(cellx, celly,...
    cellz, SingleCellDims,...
    inclusionSize, diel)
    N = size(eps_r);

    %% Set up the current source density.
    Mz = zeros(N); My = Mz; Mx = Mz;
    ind_src = [5 5 5];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
    JCurrentVector = [Mx; My; Mz];

    %% Wonsoek's scalar parameter 1, -1, or 0
    s = -1;
    %=========================================
    tic
    [A,b,Ao, bo, omega, c0, Tepsuper, Tmusuper] = ...
        solve3D_EigenEngine_Matrices(L0, wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml,s);
    toc
    %MultiCellReorder(A, b,xCells, yCells, N)
    %% parameters for a multicellreordering
    tic
    [Q, permutedIndices, hpart, vpart, partitions] = ...
    IndexPermutationCubicLattice_SLS(N,borderCoords, interiorCoords)
    toc

    %% Transform Equations
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    SymB = Q*b;



%% SOLVER STARTS HERE
    %% cellIndices lists numbers from 1 to n indentifying what object is in each cell
        %cellindices is an xCellsxyCells array
    totCells = length(partitions(2:end))
    subsystemSize = partitions(2);
    hpart = partitions(1); vpart = hpart;
    %uppermost block = boundary cell partition
    % lower blocks = interior cell partition
    Abound = SymA(1:hpart, 1:hpart);
    

    Aint = SymA(hpart+1:end, vpart+1:end);
    Aib = SymA(hpart+1:end,1:vpart);
    Abi = SymA(1:hpart, vpart+1:end);
    
    bBound = SymB(1:hpart);
    bint = SymB(hpart+1:end);
   
    App = Abound;
    Avvcell = mat2cell(Aint, subsystemSize*ones(totCells,1),...
        subsystemSize*ones(totCells,1));

    Apvcell = mat2cell(Abi, hpart, subsystemSize*ones(totCells,1));
    Avpcell = mat2cell(Aib, subsystemSize*ones(totCells,1), vpart);
    bvCell = mat2cell(bint, subsystemSize*ones(totCells,1), 1);
    bp = bBound;
    Aschur = App;
    bmod = bp;
     
    %iterate through  thAve metaatoms...linear loop (x*y), but still N^2
    % loop is linear because the Avv's we want are all on a diagonal
    SchurComplements = cell(totCells+1); %% this stores all the BD^-1C
    SchurSources = cell(totCells+1);

    parfor i = 1:totCells
        Avv = Avvcell{i,i};

        Apv = Apvcell{i};
        Avp = Avpcell{i};
        bv =bvCell{i};
        Comp = Apv*(Avv\Avp);
        SchurComplements{i} = Comp;
        SchurSources{i} = Apv*(Avv\bv);
    end
    %% Construct final Schur complement
    for i = 1:totCells
       Aschur = Aschur - SchurComplements{i};
       bmod = bmod-SchurSources{i};
    end
    
end