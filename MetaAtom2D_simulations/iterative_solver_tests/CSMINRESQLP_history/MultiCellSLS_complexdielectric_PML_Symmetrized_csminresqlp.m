%% Multicell Single Partition Analysis with real dielectric + PML
close all
clear all

% Note that the csminres-qlp solver only is useful in the case of 
% a complex dielectric with no pml as this is the only complex symmetric
% case we can formulate

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.55;  % wavelength in L0
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 80;
epsilon = 12;
Sx = SingleCellSize; Sy = SingleCellSize;

for k = 1:maxCellNumber
    N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
    N0 = N; %record the original N for comparison
    Npml = [5 5];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
    xrange = k*[-1 1];  % x boundaries in L0
    yrange = k*[-1 1];  % y boundaries in L0

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
    featureDims = [SingleCellSize/4, SingleCellSize/8, SingleCellSize/6];
    [eps_air, cellIndices] =... 
        UniformLayerCellDielectricSingleLayerSep(k, k,...
        SingleCellSize, SingleCellSize, Npml,epsilon); %% ADD coe to account for PML


    %% Set up the magnetic current source density.
    Mz = zeros(N);
    ind_src = [ceil(N/5) ceil(N/5)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;
    %Mz(75, 75) = 1;
    scale = 1e-20; %% make matrix scaling look nicer
    [A, omega,b, Sxf, Syf] = solveTE_Matrices(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
    
    %% symmetrizer
    [Pl, Pr] = SCSymmetrizer2D(Sxf, Syf, Sxf, Sxf);
    A = Pl^-1*(A*Pr^-1); 
    b = Pl\b;
    
    %% Set up the MultiCell Regime
    divx = k; divy = k;
    xCells = k; yCells = k; %% these specify the number of cells we will be dividing into
    CellDimx = N(1)/divx;
    CellDimy = N(2)/divy
    
    disp('reordering')
    tic
    [SymA, SymB, Q, permutedIndices, boundaryCells, interiorCells,pmlCell, ...
    hpart, vpart, pmlxpart, pmlypart, pmlCellDict] = ... 
        MultiCellBoundaryInteriorSLS_PML_partitioned(A, b, divx, divy, SingleCellSize, N,Npml);
    toc
    %% Perform an iterative schur complement: Aschur = sum(App - X^-1CX
    
    disp('schur complement')
    tic
    [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell, InvAvvStorage] = ...
    PrecompIterativeFourBlockSchurSingleSep_PartitionedPML(SymA, SymB, hpart, vpart, ...
    xCells, yCells, N, SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart, pmlCellDict)
    toc
    
    %% Now we need to establish formalism to extract the interior solutions
    
    %% perform iterations Analysis
    tol = 1e-14;
    maxitun = length(b);
    maxitred = length(bmod);
    
    %% solve
    M = []; %preconditioner
    tic
    [x1,flag1,iter1,Miter,QLPiter,relres1,relAres,...
          Anorm,Acond1,xnorm,Axnorm,resvec1,Aresvec]  =...
          csminresqlp(A,b,tol, length(b), M);
    toc
    tic
    [x2,flag2,iter2,Miter2,QLPiter2,relres2,relAres2,...
          Anorm2,Acond2,xnorm2,Axnorm2,resvec2,Aresvec2]  = ...
          csminresqlp(Aschur,bmod,tol,length(bmod),M);
    toc
    
    %% Visualize Field Solutions
    Hz = reshape(x1, length(x1)^.5, length(x1)^.5);
    sdf = figure;
    subplot(2,1,1)
    imagesc(abs(Hz));
    %saveas(sdf, strcat('fieldPlot with',int2str(k*k),'cells.fig'));
    
    tolSol = ...
    MultiCellSchurInteriorSolGeneral(x2, Avvcell, ...
                            Avpcell, bvCell);     
    Hzrec = reshape(Q\tolSol, length(tolSol)^.5, length(tolSol)^.5)
    subplot(2,1,2);
    imagesc(abs(Hzrec));
    f = figure()
    semilogy(resvec1/norm(b), 'linewidth', 1.5)
    hold on
    semilogy(resvec2/norm(bmod), 'linewidth', 1.5)
    %title('Convergence History for Reduced and Unreduced Solves')
    xlabel('Iteration Number')
    ylabel('log_{10}(relative residual)')
    legend('unreduced', 'reduced')
%     savefig(f, strcat('ConvHist_eps=',num2str(epsilon),'.fig'))
%     saveas(f, strcat('ConvHist_eps=',num2str(epsilon),'.png'))

    
end


