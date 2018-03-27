
%% Execute a General Reordering Function
close all
%% === ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.0;  % wavelength in L0
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 40;
xmax = 41; ymax = 42;
xmin = 2; ymin = 1;
epsilon = 1;
Sx = SingleCellSize; Sy = SingleCellSize;
    
eigenExtrems = [];
directory = strcat('D:\Nathan\Documents\StanfordYearOne',...
    '\Fan Group\FDFDDocumentation\PaperFigures\Eigenvalue_Analysis')
for k = 1:1
    
    N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
    Npml = k*[0 0];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
    xrange = k*[-1 1]*2;  % x boundaries in L0
    yrange = k*[-1 1]*2;  % y boundaries in L0

    %% Note on grid resolution of the system
    % dx/(wvlen) ~1/20 or smaller

    [xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
    %FINAL GRID PARAMETERS ARE DETERMINED AT THIS POINT
    %xrange and yrange are slightly larger
    %% Cell Division Setup
    M= N(1)*N(2);

    %% NOTE dL is not in SI units when it comes out of domain_with_pml; for our purposes, that is okay
    %for now
    resolutionFactor = max([dL(1)/N(1) dL(2)/N(2)]); %dx/N ~meters
    %spatially, what is the smallestlength scale that we have to resolve
    Nx = N(1); Ny = N(2);

    %% Set up the permittivity.
    [eps_air, cellIndices] =... 
        UniformLayerCellDielectricSingleLayerSep(k, k, ...
        SingleCellSize, SingleCellSize, Npml,epsilon); %% ADD coe to account for PML
    
    %% put inclusion inside schur
    %eps_air(5:15,10:30) = 4;
    %% put it outside
    %eps_air(25:35,15:35) = 4i;
    %% Set up the magnetic current source density.
    Mz = zeros(N);fi
    ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;
    scale = 1e-14; %% make matrix scaling look nicer
    [A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices_unscaled(wvlen, xrange, yrange, eps_air, Mz, Npml);
    A = scale*A; b = scale*b;
    
    %% Once the DOF Are Selected, we execute the reordering function
    fig = figure;
    imagesc(abs(eps_air));
%     p = ginput(2)

    %% Do the reordering of the selected coordinates
    % the input will be two points, which is enough to specify a rectangle
    
%     sp(1) = min(floor(p(1)), floor(p(2))); %xmin
%     sp(2) = min(floor(p(3)), floor(p(4))); %ymin
%     sp(3) = max(ceil(p(1)), ceil(p(2)));   %xmax
%     sp(4) = max(ceil(p(3)), ceil(p(4)));   %ymax
    
    [NaturalOrdering, OrderMap] = CoordinateIndexing(N(1), N(2));
    
    %% artifical sp
    sp = [xmin,ymin, xmax, ymax];
    
    %% visualize Schur complement
    map = zeros(Nx, Nx);
    map(sp(1):sp(3), sp(2):sp(4)) = 1;
    map(1:sp(1)-1,:) = 2;
    map(sp(3)+1:end, :) = 3;
    map(sp(1):sp(3), 1:sp(2)) = 4;
    map(sp(1):sp(3), sp(4):end) = 5;
    figure;
    imagesc(map)

    %% Create Storage Data Structures
    interior = OrderMap(sp(1):sp(3), sp(2):sp(4));
    
    exterior1 = OrderMap(1:sp(1)-1,:);
    exterior3 = OrderMap(sp(3)+1:end, :);
    exterior2 = OrderMap(sp(1):sp(3), 1:sp(2));
    exterior4 = OrderMap(sp(1):sp(3), sp(4):end);
    
    %% reshape exteriors
    exterior1 = reshape(exterior1, numel(exterior1),1);
    exterior2 = reshape(exterior2, numel(exterior2),1);
    exterior3 = reshape(exterior3, numel(exterior3),1);
    exterior4 = reshape(exterior4, numel(exterior4),1);
    
    interiorCellIndex = reshape(interior, numel(interior),1);
    exteriorCellIndex = [exterior1;exterior2; exterior3; exterior4];
    permutedIndices = [exteriorCellIndex; interiorCellIndex;];
    
    %% Partition Sizes
    hpart = length(exteriorCellIndex);
    vpart = hpart;
    
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
    SymB = Q*b;
        
    %% Test the Schur complement (Fourblockschur)
    [Aschur, bmod, App, Apv, Avp, Avv, bp, bv] = ...
    FourblockSchur(SymA,SymB, hpart, vpart);
    
    %% permute Aschur (since it is sparse)
    r = symrcm(Aschur);
    Aschur_b = Aschur(r,r);
    m = symamd(Aschur);
    Aschur_m = Aschur(m,m);
    
    
    %% verify solutions
    figure;
    plot(abs(Q.'*(SymA\SymB))); 
    hold on
    plot(abs(A\b));
    title('Permutation Soln Check')
    figure;
    plot(abs(SymA\SymB));
    hold on;
    plot(abs(Aschur\bmod));

    %% test Aschur modes;
    for m = 1:5
        [Vs, Ds, Ws] = eigs(Aschur, m, 'sm');
%         v = reshape(Vs(:,m),xmax-xmin+1, ymax-ymin+1);
%         figure;
%         visreal(v, xrange, yrange);
    end
    
    %% test full modes
    [Vf, Df, Wf] = eigs(A, 1, 'sm');
    vf = reshape(Vf, length(eps_air), length(eps_air));
    figure;
    visreal(vf, xrange, yrange)
    
    [V1, D1, W1] = eig(full(Aschur));
    [V2, D2, W2] = eig(full(A));
    
    figure;
    semilogy(sort(abs(diag(D1)))); 
    hold on;
    semilogy(sort(abs(diag(D2))));

    %% quick iterations test
    [xf, flagf, relresf, iterf, resvecf] = qmr(A, b, 1e-10, length(b));
    [xs, flags, relress, iters, resvecs] = qmr(Aschur, bmod, 1e-10, length(bmod));
    figure;
    semilogy(resvecf);
    hold on;
    semilogy(resvecs);
    
    %% Analyze the App block...is a dirichlet boundary condition operator?
    bs = zeros(31^2, 1); bs(480) = 1i;
    soln = App\bs;
    sol = reshape(soln, 31, 31);
    
    [Vpp, Dpp, Wpp] = eig(full(App));
    figure;
    visreal(reshape(Vpp(:,end), 31, 31), xrange, yrange);
    
    [Vvv, Dvv, Wvv] = eig(full(Avv));
    figure;
    visreal(reshape(Vvv(:,end), 11, 73), xrange, yrange);
    
end