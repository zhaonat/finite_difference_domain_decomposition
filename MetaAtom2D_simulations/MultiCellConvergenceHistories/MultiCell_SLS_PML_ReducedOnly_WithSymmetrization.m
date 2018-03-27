%% Multicell Single Partition Analysis with real dielectric + PML
close all
clear 

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.55*L0;  % wavelength in L0
maxCellNumber = 1;
SingleCellSize = 70;
epsilon = 12;
Sx = SingleCellSize; Sy = SingleCellSize;

%% ==================== Create Data Directory
directory = strcat('D:\Nathan\Documents\StanfordYearOne\Fan Group\FDFDSchurComplement\MetaAtom2D\',...
                    'MultiCellsSingleLayerSeparator_PartitionedPML');

dims = [];
for k = 3:3
    for epsilon = [12i]
        g = k;
        N = [g*SingleCellSize+g+1 k*SingleCellSize+k+1];  % [Nx Ny]
        N0 = N; %record the original N for comparison
        Npml = [25,25];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
        xrange = g*L0*[-1 1];  % x boundaries in L0
        yrange = k*L0*[-1 1];  % y boundaries in L0

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
        ind_src = [Npml(1)+SingleCellSize Npml(2)+SingleCellSize];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
        Mz(ind_src(1), ind_src(2)) = 1;
        %Mz(75, 75) = 1;
        %scale = 1e-13; %% make matrix scaling look nicer
        [A, omega,b, Sxf, Syf, Sxb, Syb] = solveTE_Matrices_unscaled(wvlen, xrange, yrange, eps_air, Mz, Npml);
        %A = scale*A; b = scale*b;
        
        %% symmetrize the system
        [Pl, Pr] = SCSymmetrizer2D(Sxf, Syf, Sxb, Syb)
        A0 = A; b0 = b;
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
        
        %% Do it with the unsymmetrized system
                disp('reordering')
        tic
        [SymA0, SymB0, Q0, permutedIndices0, boundaryCells0, interiorCells0,pmlCell0, ...
        hpart0, vpart0, pmlxpart0, pmlypart0, pmlCellDict0] = ... 
            MultiCellBoundaryInteriorSLS_PML_partitioned(A0, b0, divx, divy, SingleCellSize, N,Npml);
        toc
        %% Perform an iterative schur complement: Aschur = sum(App - X^-1CX

        disp('schur complement')
        tic
        [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvcell, InvAvvStorage] = ...
            PrecompIterativeFourBlockSchurSingleSep_PartitionedPML(SymA, SymB, hpart, vpart,divx, divy,N,...
            SingleCellSize, cellIndices, Npml, pmlxpart, pmlypart, pmlCellDict);
        toc
        
        disp('schur complement')
        tic
        [Aschur0, bmod0] = ...
            PrecompIterativeFourBlockSchurSingleSep_PartitionedPML(SymA0, SymB0, hpart, vpart,divx, divy,N,...
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
%         disp('full solve');
%         tic
%         [xTol, flag2, relres2, iter2, resvec2, iterates2] = ...
%             qmr_iterates(A, b, tol, maxitun, [],[],[],Sampling2);
%         toc
        
        %% now we recover the residuals on the iterates
        % recover the tol sol
%         disp('iterate recovery')
%         s1 = size(iterates1); %second dimension is the number of iterates
%         relresRedTot = [];
%         for i = 1:s1(2)
%             i
%             schur_it = iterates1(:, i);
%             tolSol = ...
%             MultiCellSchurInteriorSolGeneral(schur_it, Avvcell, ...
%                             Avpcell, bvcell);  
%             relresRedTot = [relresRedTot, norm(SymA*tolSol - SymB)];
%         end
        
        %% reconstruct the tolSol
        tolSol = ...
            MultiCellSchurInteriorSolGeneral(xSchur, Avvcell,  ...
                            Avpcell, bvcell);     
        figure; 
        Hzrec = reshape(Q\tolSol, length(tolSol)^.5, length(tolSol)^.5);
        visreal(Hzrec, xrange, yrange);
        save(strcat('Npml=',num2str(Npml(1)),'_numcells=',num2str(k),'_eps=',num2str(epsilon),...
            '_MultiCellPMLConvergenceHistories.mat'));
       % 'iterations1', 'iterations2',
       %     'resvec1', 'resvec2', 'bmod', 'b', 'relresBound', 'relresRedTot');

        xdirect = A\b;
        Hz = reshape(xdirect, Nx, Ny);
        visreal(Hz, xrange, yrange);
    end
    
end






