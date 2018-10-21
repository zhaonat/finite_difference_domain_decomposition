%% Multicell Single Partition Analysis with real dielectric + PML
close all
clear all

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 3.0;  % wavelength in L0
iterScaling = []; 
maxCellNumber = 1;
SingleCellSize = 80;
epsilon = -9.8+0.31i;
Sx = SingleCellSize; Sy = SingleCellSize;

%% ==================== Create Data Directory

for npml = [15, 0]
    for k = 8:9
        N = [k*SingleCellSize+k+1 k*SingleCellSize+k+1];  % [Nx Ny]
        N0 = N; %record the original N for comparison
        Npml = [npml, npml];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
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
        featureDims = SingleCellSize*[1/4, 1/4, 1/6];
        [eps_air, cellIndices] =... 
            multiRandomCellDielectricSingleLayerSep(k, k, SingleCellSize, ...
            SingleCellSize, Npml,epsilon, featureDims); %% ADD coe to account for PML


        %% Set up the magnetic current source density.
        Mz = zeros(N);
        ind_src = [ceil(N/2) ceil(N/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
        Mz(ind_src(1), ind_src(2)) = 1;
        %Mz(75, 75) = 1;
        %scale = 1e-13; %% make matrix scaling look nicer
        [A, omega,b, Sxf, Dxf,Dyf] = solveTE_Matrices_unscaled(wvlen, xrange, yrange, eps_air, Mz, Npml);
        %A = scale*A; b = scale*b;

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
        [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvcell, InvAvvStorage] = ...
            PrecompIterativeFourBlockSchurSingleSep_PartitionedPML(SymA, SymB, hpart, vpart,divx, divy,N,...
            SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart, pmlCellDict);
        toc

        %% Now we need to establish formalism to extract the interior solutions

        %% perform iterations Analysis
        tol = 1e-14;
        maxitun = length(b);
        maxitred = length(bmod);
        %% We need to extract the convergence history for the border and the interior
        disp('full solve')
        tic 
        [x1, flag1, relres1, iter1, resvec1] = qmr(SymA, SymB, tol, maxitun);
        toc
        relres1
        iter1
        disp('reduced solve')
        tic
        [x2, flag2, relres2, iter2, resvec2] = qmr(Aschur, bmod, tol, maxitred);
        toc
        iterScaling = [iterScaling;iter1, iter2]
        %% Visualize Field Solutions
        Hz = reshape(Q\x1, length(x1)^.5, length(x1)^.5);
        f = figure;
        subplot(2,2,1)
        visreal(Hz, xrange, yrange);
        title('Unreduced Field Solution')
        hcb=colorbar;
        high = full(max(max(real(Hz))));
        low = full(min(min(real(Hz))));
        set(hcb, 'YTick', [low,0, high])
        set(hcb,'TickLabels',{'-max', '0', 'max'})
        set(gca, 'LineWidth',1.25)

        tolSol = ...
        MultiCellSchurInteriorSolGeneral(x2, Avvcell, ...
                                Avpcell, bvcell);     
        Hzrec = reshape(Q\tolSol, length(tolSol)^.5, length(tolSol)^.5);
        subplot(2,2,3);
        visreal((Hzrec), xrange, yrange);
        title('Schur Field Solution')
        hcb=colorbar;
        high = full(max(max(real(Hz))));
        low = full(min(min(real(Hz))));
        set(hcb, 'YTick', [low,0, high])
        set(hcb,'TickLabels',{'-max', '0', 'max'})
        set(gca, 'LineWidth',1.25)


        subplot(2,2,[2,4])
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
    fileName = strcat('Multicell2D_wvlen=',num2str(3),'Npml=',num2str(npml),...
        '_epsilon=',num2str(epsilon),'_ConvergenceHistory.mat');
    save(fileName,'resvec1', 'resvec2', 'b', 'bmod', 'flag1', 'flag2', 'relres1', 'relres2',...
        'iter1', 'iter2', 'x1', 'x2');
    
end
% tic
% t1 = Aschur\bmod;
% toc
% tic
% t2 = A\b;
% toc

%% Plot cell scaling properties
% g = figure;
% plot(iterScaling)
% title('Cell Scaling for PML (fully complemented out)')
% xlabel('Iteration #')
% ylabel('log_{10}(Relative R1esidual)')
% saveas(g, 'IterationScalingWithPML (fully complemented out).fig')



