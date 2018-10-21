close all
clear 

%% Demonstration of General Domain Partitioning
L0 = 1e-6;  % length unit: microns
wvlen = 1;  % wavelength in units of L0, DO NOT MULTIPLY!
epsilon = 12;

N = [100,100];  % [Nx Ny]
N0 = N; %record the original N for comparison
Npml = [15 15];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
xrange = [-1 1];  % x boundaries in L0
yrange = [-1 1];  % y boundaries in L0

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PM
%% Cell Division Setup
M= N(1)*N(2);

%% epsilon
epsilon = ones(N);

%% source
Mz = zeros(N);
Mz(50,50) = 1;

%% Get system matrices
[A, b, omega] = solveTE_matrices(L0, wvlen, xrange, ...
    yrange, epsilon, Mz, Npml);

%% add pec
PEC_mask_x = create_PEC(N, 'x');
PEC_mask_y = create_PEC(N, 'y');
%A = PEC_mask_x*PEC_mask_y*A*PEC_mask_y*PEC_mask_x;
%% define a subdomain rectangle
coord1 = [40,40];
coord2 = [90,90];

exist_mask = ones(N);
[Q,mask] = square_partition(N, coord1, coord2, exist_mask);

%% visualize mask
figure();
spy(mask)

%% permute system
SymA = Q*A*Q.';
SymB = Q*b;

%% test solutions
tic
x0 = A\b;
toc
x = Q\(SymA\SymB);
figure(); plot(real(x)); hold on; plot(real(x0))
hold on;
%% for kicks, demonstrate the four block schur function
hpart =  M-nnz(mask);
tic
[Aschur, bmod, App, Apv, Avp, Avv, bp, bv] = ...
    FourblockSchur(SymA,SymB, hpart);
toc
tic
xSchur = Aschur\bmod;
toc

%% reconstruct
recField = reconstruct_four_block(xSchur, Avv, Avp, bv, Q);

Hz = reshape(recField, N(1), N(2));
figure();visreal(Hz, xrange, yrange);

%% in the end, we will need to figure out how to calculate the adjoint
