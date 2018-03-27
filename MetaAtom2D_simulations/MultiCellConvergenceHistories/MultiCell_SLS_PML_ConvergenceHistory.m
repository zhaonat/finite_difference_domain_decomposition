%% Multicell Single Partition Analysis with real dielectric + PML
close all
clear 

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 3;  % wavelength in L0
SingleCellSize = 80;
Sx = SingleCellSize; Sy = SingleCellSize;

dims = [];
for wvlen = [3];
for k = 8:8
    for epsilon = [12i]
        g = k;
        N = [g*SingleCellSize+g+1 k*SingleCellSize+k+1];  % [Nx Ny]
        N0 = N; %record the original N for comparison
        Npml = [15,15];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
        xrange = g*[-1 1];  % x boundaries in L0
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
        featureDims = [SingleCellSize/2, SingleCellSize/4, SingleCellSize/6];
        [eps_air, cellIndices] =... 
            multiRandomCellDielectricSingleLayerSep(k, g, SingleCellSize,...
            SingleCellSize, Npml,epsilon, featureDims); %% ADD coe to account for PML
        

        %% Set up the magnetic current source density.
        Mz = zeros(N);
        ind_src = [round(Nx/2), round(Ny/2)];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
        Mz(ind_src(1), ind_src(2)) = 1;
        %Mz(75, 75) = 1;
        %scale = 1e-13; %% make matrix scaling look nicer
        [A, omega,b, Sxf, Syf, Sxb, Syb] = solveTE_Matrices(L0,wvlen, xrange, yrange, eps_air, Mz, Npml);
        %A = scale*A; b = scale*b;
        
        %% symmetrize the system
        [Pl, Pr] = SCSymmetrizer2D(Sxf, Syf, Sxb, Syb);
        A = Pl^-1*(A*Pr^-1);
        b = Pl\b;
        
        %% Set up the MultiCell Regime
        divx = g; divy = k;
        xCells = g; yCells = k; %% these specify the number of cells we will be dividing into

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
        
        %% Convergence Histories
        %% QMR
        Sampling1 = round(length(Aschur)/100);
        Sampling2 = round(length(A)/600);
        disp('reduced solve')
        tic 
        [xSchur, flag1, relres1, iter1, resvec1, iterates1] = ...
            qmr_iterates(Aschur, bmod, tol, maxitred, [],[],[],Sampling1);
        toc
        figure; semilogy(resvec1);
        disp('full solve');
        tic
        [xTol, flag2, relres2, iter2, resvec2, iterates2] = ...
            qmr_iterates(A, b, tol, maxitun, [],[],[],Sampling2);
        toc
        
        %% now we recover the residuals on the iterates
        % recover the tol sol
        disp('iterate recovery')
        s1 = size(iterates1); %second dimension is the number of iterates
        relresRedTot = [];
        for i = 1:s1(2)
            i
            schur_it = iterates1(:, i);
            tolSol = ...
            MultiCellSchurInteriorSolGeneral(schur_it, Avvcell, ...
                            Avpcell, bvcell);  
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
        semilogy(iterations1, relresRedTot/norm(b), 'linewidth', 1.5)
        semilogy(iterations2, relresBound/norm(bmod), 'linewidth', 1.5)
        semilogy(resvec2/norm(b), 'linewidth', 1.5)
        %title('Convergence History for Reduced and Unreduced Solves')
        xlabel('Iteration Number')
        ylabel('log_{10}(relative residual)')
        legend('r1', 'r2', 'ur1', 'ur2')
        
        %% reconstruct the tolSol
        tolSol = ...
            MultiCellSchurInteriorSolGeneral(xSchur, Avvcell, ...
                            Avpcell, bvcell);     
        figure; 
        Hzrec = reshape(Q\tolSol, length(tolSol)^.5, length(tolSol)^.5);
        visreal(Hzrec, xrange, yrange);
        save(strcat('wvlen=',num2str(wvlen),'_Npml=',num2str(Npml(1)),'_numcells=',num2str(k),'_eps=',num2str(epsilon),...
            '_MultiCellPMLConvergenceHistories.mat'));
       % 'iterations1', 'iterations2',
       %     'resvec1', 'resvec2', 'bmod', 'b', 'relresBound', 'relresRedTot');


    end
    
end
end








