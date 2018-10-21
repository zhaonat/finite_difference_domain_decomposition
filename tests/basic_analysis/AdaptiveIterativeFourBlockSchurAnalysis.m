%% MAIN IDEA, BE ABLE TO FIND THE CORRECT INTERIOR PARTITIONS BY THE MATRIX STRUCTURE ALONE
clear; close all; clc;

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 10.0*L0;  % wavelength in L0
xrange = L0*[-5 5];  % x boundaries in L0
yrange = L0*[-5 5];  % y boundaries in L0
N = [40 40];  % [Nx Ny]
M= N(1)*N(2);
Npml = [0 0];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0

%% Note on grid resolution of the system
% dx/(wvlen) ~1/20 or smaller

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
%FINAL GRID PARAMETERS ARE DETERMINED AT THIS POINT
%xrange and yrange are slightly larger

%% NOTE dL is not in SI units when it comes out of domain_with_pml; for our purposes, that is okay
%for now
resolutionFactor = max([dL(1)/N(1) dL(2)/N(2)]); %dx/N ~meters
%spatially, what is the smallestlength scale that we have to resolve

%% Set up the permittivity.
eps_air = ones(N);
Nx = N(1); Ny = N(2);
for i = 1:N
   for j = 1:N
       if(j >10 && j<20)
          eps_air(i,j)= 1; 
       end
   end
end

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;
%Mz(75, 75) = 1;

%% Solve TE equations.
[Hz, Ex, Ey, A, omega,b, Sxf, Dxf,Dyf, sxf, syf] = solveTE_unscaled(wvlen, xrange, yrange, eps_air, Mz, Npml);
figure;
visabs(Hz, xrange, yrange);

xCells = 2; yCells = 2;
%% Execute Symmetry Operation
% [Q, SymA, SymB, permutedIndices, boundaryCells, interiorCells, hpart, vpart] = ...
%     MultiCellBoundaryInteriorSinglePartition(A, b, xCells, yCells, N);
[SymA, SymB, Q, permutedIndices, boundaryInd, interiorInd, hpart, vpart] = ...
    MultiCellReorder(A,b, xCells, yCells, N);

%% PRESENTLY THERE IS A PROBLEM IN the AVVCells, the interiors
%% aren't fully decoupled from each other....

%% Solver Code Starts Here
CDx = Nx/xCells; CDy = Ny/yCells;
totCells = xCells*yCells;


    %uppermost block = boundary cell partition
    % lower blocks = interior cell partition
    Abound = SymA(1:hpart, 1:vpart);
    

    Aint = SymA(hpart+1:end, vpart+1:end);
    Aib = SymA(hpart+1:end,1:vpart);
    Abi = SymA(1:hpart, vpart+1:end);
    
    bBound = SymB(1:hpart);
    bint = SymB(hpart+1:end);
    Nx = N(1); Ny = N(2);
    App = Abound;
    Avv = cell(1,1);
    %% now we have to find the breaks in Aint
    start = 1;
    count = 1;
    for i = 1:length(Aint)
       if(nnz(Aint(:,i))==)
           endind = i;
           Avv{count, count} = Aint(start:endind, start:endind);
           start = endind;
           count = count+1;
       end
    end
    
    %% this formulation is FASTER...
    
%% Compare to Single Block Partition
to = cputime;
[AschurWhole, bmodWhole] = FourblockSchur(SymA, SymB, hpart, vpart);
trun2 = cputime-to
sol2 = AschurWhole\bmodWhole;

figure
plot(abs(sol2))
figure;
plot(abs(Aschur\bmod))