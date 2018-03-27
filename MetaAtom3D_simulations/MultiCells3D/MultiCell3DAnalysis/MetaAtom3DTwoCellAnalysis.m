clear;
close all;
cd('D:\Nathan\Documents\StanfordYearOne\Fan Group\FDFDSchurComplement\MetaAtom3D\MultiCells3D\MultiCell3DAnalysis')

%% Parameter Set-up
L0 = 1e-6;  % length unit: microns
wvlen = 30.0;  % wavelength in L0
xSize = 2;
xrange = xSize*[-5 5];  % x boundaries in L0
yrange = [-5 5];  % y boundaries in L0
zrange = [-5 5];
CellDim = 10;
Nx = xSize*CellDim; Ny = CellDim; Nz = CellDim;
N = [Nx Ny Nz];  % [Nx Ny]
M = N(1)*N(2)*N(3);
Npml = [0 0 0];  % [Nx_pml Ny_pml]

%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
eps_r = ones(N);
eps_r(4:7, 4:7, 4:7) = 20; eps_r(12:18, 2:8, 2:8) = 12;
%% Set up the current source density.
Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = [5 5 5];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
JCurrentVector = [Mx; My; Mz];
%% Wonsoek's scalar parameter 1, -1, or 0
s = -1;
%=========================================

[A, b, Ao, bo, omega, c0, TepsSuper, TmuSuper]= ...
    solve3D_EigenEngine_Matrices(L0, wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml,s);

tic
[permutedIndices, boundaryCells, interiorCells, hpart, vpart, l1, l2] = ...
    IndexPermutation3DTwoCell(N,2);
toc
cellpartitions = [l1, l2];

%% visualize the partition
figure;
g1 = length(boundaryCells);
plot3(boundaryCells(:,1), boundaryCells(:,2), boundaryCells(:,3), '.g', 'markersize', 10);
hold on;
g = length(interiorCells);
plot3(interiorCells(1:l1,1), interiorCells(1:l1,2), interiorCells(1:l1,3), '.b', 'markersize', 5)
hold on;
plot3(interiorCells(l1+1:end,1), interiorCells(l1+1:end,2), interiorCells(l1+1:end,3), '.r', 'markersize', 5)

xind = zeros(3*Nx*Ny*Nz,1);
yind = zeros(3*Nx*Ny*Nz,1);
vals = ones(3*Nx*Ny*Nz,1);
for i = 1:3*N(1)*N(2)*N(3)

   indexshift = permutedIndices(i);
   xind(i) = i;
   yind(i) = indexshift;
   %Q(i,indexshift) = 1;
end
Q = sparse(xind,yind,vals);

SymA = Q*A*transpose(Q);
SymB = Q*b;

% x2 = SymA\symB;
NaturalOrdering = CoordinateIndexing3D(N(1), N(2), N(3));
% x1 = A\b; % solutions are in fact the same, which is very good

%% schur complement is now positive definite with Wonsoek's accelerator
%% Execute Schur Complement
part1 = hpart;
%% adaptive iterative four block schur
tic
[Aschur, bmod, App, Apv, Avp, Avv, bp, bv] = ...
    FourblockSchur(SymA, SymB, part1, part1);
toc
%% reconstruct interior solution (is sol2)
%recSol = [sol1; sol2];

%% do a quick iterative test
% tic
% [xTol, flag1, relres1, iter1, resvec1] = minres(A,b, 1e-14, 1e6);
% toc
% tic
% [xSchur, flag2, relres2, iter2,resvec2] = minres(Aschur,bmod, 1e-14, 1e6);
% toc
tic
xTol = Ao\bo; %%without wonsoek's accelerator, so nnz(matrix) much higher
toc
tic
xSchur = Aschur\bmod;
toc
%% Visualize Matrices
As = figure;
spy(Aschur)
saveas(As, strcat('Aschur for cellDim = ', num2str(Ny)))

Bs = figure;
spy(SymA)
saveas(Bs, strcat('SymA for cellDim = ', num2str(Ny)))


% %% visualize resvecs
% figure;
% plot(log10(resvec1/norm(b)));
% hold on;
% plot(log10(resvec2/norm(bmod)));
%% Visualize Field
%[Exc, Eyc, Ezc] = FDFD3D_SliceVisualization(Ex,Ey, Ez, N, xrange, yrange);

%% Verify solutions on the boundary are consistent with the original solution
% 
% ExOriginalBoundary = VectorBoundaryPoint3DExtraction(Ex,N);
% EyOriginalBoundary = VectorBoundaryPoint3DExtraction(Ey,N);
% EzOriginalBoundary = VectorBoundaryPoint3DExtraction(Ez,N);
% boundDim = length(Ex);
% ExSchurBound = sol1(1:400); %%this has to be rearranged to get the boundary
% figure;
% plot(abs(ExOriginalBoundary));
% title('original boundary field')
% figure;
% plot(abs(ExSchurBound));
% title('extracted schur solution boundary')
% figure;
% spy(Aschur)
% figure;
% imagesc(log(abs(Aschur)))

%% Do some iterations analysis
tol = 1e-11;
Sampling = 2;
maxit = 1e6;
[xTol, flag, iter1, relresTot, relresBound, xBound] =...
NZqmrUnReduced( SymA, SymB, tol, maxit, M, Aschur, bmod,hpart, Sampling)

[xSchur, flag, iter, relresSchur, relresRedTot] =...
NZqmrReducedSingleCell( Aschur, bmod, tol, maxit,  SymA, SymB, Avv, Avp, Sampling, [], hpart)
%% sum( sum(abs(Aschur - Aschur.')))
%plot(abs(sol1));
%conditionNumbers = [conditionNumbers; cond(A) cond(Aschur)];
%     iterate0_array = [iterate0_array, iter1]; 
%     iterateSchur_array = [iterateSchur_array iter2];
%     
convhist = figure;
plot(log10(relresTot), 'bo');
hold on;
plot(log10(relresBound), 'rx');
plot(log10(relresSchur), 'r.');
plot(log10(relresRedTot), 'b.')
legend('total residual', 'residual on boundary', 'residual schur bound', 'residual total schur')
xlabel('iterations')
ylabel('log10(relres)')
title('Iterations Analysis for 3D Single Cells')
saveas(convhist, strcat('Convergence History for Nx=',num2str(i)));