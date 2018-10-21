%% Single layer Schur complement with PML setup
close all
clear

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.55*L0;  % wavelength in L0
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 40;
epsilon = 1;
Sx = SingleCellSize; Sy = SingleCellSize;

k = 2

%% ================ Domain Size Parameters ==========================
N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
N0 = N;
Npml = [10 10];  
xrange = k*L0*[-1 1];  % x boundaries in L0
yrange = k*L0*[-1 1];  % y boundaries in L0

%% Generate Parameters for Domain with PML
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  

%% Generate Parameters for same domain with No PML
[xrange1, yrange1, N1, dL1, Lpml1] = domain_with_pml(xrange, yrange, N, [0,0]);  


%% ==============Setup Dielectric ===================================
[eps_air, cellIndices] = ... 
    multiRandomCellDielectricSingleLayerSep(k, k, SingleCellSize,...
    SingleCellSize, Npml,epsilon); %% ADD coe to account for PML

[eps_air1, cellIndices1] = ... 
    multiRandomCellDielectricSingleLayerSep(k, k, SingleCellSize,...
    SingleCellSize, [0 0],epsilon); %% ADD coe to account for PML

imagesc(abs(eps_air))

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

%no pml, matrices have the same size...
[A1, omega1,b1, Sxf1, Dxf1,Dyf1] = solveTE_Matrices_unscaled(wvlen, xrange1, ...
    yrange1, eps_air1, Mz1, [0 0]);

%% Test Solve
x = A\b;
Hz = reshape(x, length(x)^.5, length(x)^.5);
visabs(abs(Hz), xrange, yrange)


