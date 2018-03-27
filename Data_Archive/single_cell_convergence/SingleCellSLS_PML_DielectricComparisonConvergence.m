%%Single Cell Analysis with vacuum, real dielectric, and complex dielectric + PML
close all
clear
%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 80;
Sx = SingleCellSize; Sy = SingleCellSize;

k = 1;

counter = 1
for Nl = [0]
for wvlen = [0.1]
for epsilon = [1,12,12i]
    solns = cell(2,1);
    resvecs = cell(4,1);
    N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
    N0 = N; %record the original N for comparison
    Npml = [Nl Nl];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
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
    [eps_air, cellIndices] =... 
        multiRandomCellDielectricSingleLayerSep(k, k, SingleCellSize, SingleCellSize, Npml,epsilon); %% ADD coe to account for PML


    %% Set up the magnetic current source density.
    Mz = zeros(N);
    ind_src = [ceil(N/5) ceil(N/5)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;
    %Mz(75, 75) = 1;
    %scale = 1e-13; %% make matrix scaling look nicer
    [A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
    %A = scale*A; b = scale*b;
    
    %% Set up the MultiCell Regime
    divx = k; divy = k;
    xCells = k; yCells = k; %% these specify the number of cells we will be dividing into
    
    disp('reordering')
    tic
    [SymA, SymB, Q, permutedIndices, boundaryCells, interiorCells,pmlCell, ...
    hpart, vpart, pmlxpart, pmlypart] = ... 
        MultiCellBoundaryInteriorSLS_PML(A, b, divx, divy, SingleCellSize, N,Npml);
    toc
    %% Perform an iterative schur complement: Aschur = sum(App - X^-1CX
    
    disp('schur complement')
    tic
    [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvcell, InvAvvStorage] = ...
        PrecompIterativeFourBlockSchurSingleSep_PML(SymA, SymB, hpart, vpart,divx, divy,N,...
        SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart);
    toc
    

    %% perform iterations Analysis

    %% We need to extract the convergence history for the border and the interior
    %% GMRES

    tol = 1e-14;
    %% QMR
    Sampling1 = 10;
    Sampling2 = 20;
    relresRedtot = []; relresBound = [];
    disp('reduced solve')
    tic 
    [xSchur, flag1, relres1, iter1, resvec1, iterates1] = ...
        qmr_iterates(Aschur, bmod, tol, length(bmod), [],[],[],Sampling1);
    toc
    disp('full solve');
    tic
    [xTol, flag2, relres2, iter2, resvec2, iterates2] = ...
        qmr_iterates(A, b, tol, length(b), [],[],[],Sampling2);
    toc

    %% now we recover the residuals on the iterates
    % recover the tol sol
    disp('iterate recovery')
    s1 = size(iterates1); %second dimension is the number of iterates
    relresRedTot = [];
    [L,U] = lu(Avvcell{1,1});
    bv = bvcell{1};
    for i = 1:s1(2)
        schur_it = iterates1(:, i);
        tolSol = [schur_it; U\(L\(bv - Avpcell{1}*schur_it))];
        relresRedTot = [relresRedTot, norm(SymA*tolSol - SymB)];
    end

    s2 = size(iterates2); %second dimension is the number of iterates
    relresBound = [];
    for i = 1:s2(2)
        tol_it = Q*iterates2(:, i); %solution is done with A --> to go SymA
        bound = tol_it(1:hpart); 
        relresBound = [relresBound, norm(Aschur*bound - bmod)];
    end

    %% semilogy
    iterations1 = linspace(1, length(resvec1), length(relresRedTot));
    iterations2 = linspace(1, length(resvec2), length(relresBound));
    a = figure;
    semilogy(resvec1/norm(bmod), 'linewidth', 1.5)
    hold on
    semilogy(iterations1(1:end), relresRedTot/norm(b), 'linewidth', 1.5)
    semilogy(iterations2(1:end), relresBound/norm(bmod), 'linewidth', 1.5)
    semilogy(resvec2/norm(b), 'linewidth', 1.5)
    %title('Convergence History for Reduced and Unreduced Solves')
    xlabel('Iteration Number')
    ylabel('log_{10}(relative residual)')
    legend('r1', 'r2', 'ur1', 'ur2')
    flag1
    flag2
    %% Visualize Field Solutions
    Hz = reshape(Q*xTol, length(xTol)^.5, length(xTol)^.5);
    sdf = figure;
    subplot(2,2,1)
    cmax = visreal(i*Hz, xrange, yrange);
    title('Unreduced Field Solution')
    hcb=colorbar;
    
    set(hcb, 'YTick', [-cmax,0, cmax])
    set(hcb,'TickLabels',{'-max', '0', 'max'})
    set(gca, 'LineWidth',1.25)

    tolSol = ...
    MultiCellSchurInteriorSolGeneral(xSchur, Avvcell, InvAvvStorage, ...
                            Avpcell, bvcell);     
    Hzrec = reshape(Q\tolSol, length(tolSol)^.5, length(tolSol)^.5)
    
    %% STORAGE
    solns{ 1} = Hzrec;
    solns{ 2} = Hz;

    resvecs{1,2} = resvec1;
    resvecs{2,2}= resvec2;
    resvecs{3,2}= relresRedTot;
    resvecs{4,2}= relresBound;
    resvecs{1,1} = 1:length(resvec1);
    resvecs{2,1}= 1:length(resvec2);
    resvecs{3,1}= iterations1;
    resvecs{4,1}= iterations2;
    
    %% VISUALIZATION
    subplot(2,2,3);
    visreal(i*(Hzrec), xrange, yrange);
    
    title('Schur Field Solution')
    hcb=colorbar;
    high = full(max(max(real(Hz))));
    low = full(min(min(real(Hz))));
    %set(hcb, 'YTick', [low,0, high])
    %set(hcb,'TickLabels',{'-max', '0', 'max'})
    set(gca, 'LineWidth',1.25)
    
    subplot(2,2,[2,4])
    semilogy(resvec1/norm(bmod), 'linewidth', 1.5)
    hold on
    semilogy(iterations1, relresRedTot/norm(b), 'linewidth', 1.5)
    semilogy(iterations2, relresBound/norm(bmod), 'linewidth', 1.5)
    semilogy(resvec2/norm(b), 'linewidth', 1.5)
    %title('Convergence History for Reduced and Unreduced Solves')
    xlabel('Iteration Number')
    ylabel('log_{10}(relative residual)')
    legend('r1', 'r2', 'ur1', 'ur2')
    
    savefig(sdf, strcat('Nx=',num2str(N(1)),'wvlen=',num2str(wvlen),'_Soln+TotalConvHist_eps=',num2str(epsilon),'.fig'))
    saveas(sdf, strcat('Nx=',num2str(N(1)),'wvlen=',num2str(wvlen),'_Soln+TotalConvHist_eps=',num2str(epsilon),'.png'))
    save(strcat('Npml=',num2str(Npml(1)),'_Nx=',num2str(N(1)),'wvlen=',...
        num2str(wvlen),'_Soln+TotalConvHist_eps=',num2str(epsilon),'.mat'),...
        'solns', 'resvecs', 'b', 'bmod', 'A', 'Aschur', 'Q', 'eps_air');
    
    
end
end
end
