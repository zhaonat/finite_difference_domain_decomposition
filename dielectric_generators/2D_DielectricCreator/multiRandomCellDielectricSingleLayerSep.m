%% Dielectric MultiCell Formatter
%% a function which generates a 2D eps_dielectric matrix given the number of cells
%% in the domain and the protoypical cell

function [finalEps, cellIndices] = ...
    multiRandomCellDielectricSingleLayerSep(xCells, yCells,celldimx, ...
    celldimy, Npml, dielConst, featureDims)
    rng(1) %% set seed for repeatability so same set is generated
    %% Parameter Descriptions
    % xCells - # of metaatoms in x direction
    % yCells - # of metaatoms in y direction
    % cellindices - a dictionary of what dielectric inclusions is
    %     contained in which meatatom; this is a matrix with the dimension of
    %     the metatom grid: xcells*ycells;
    % finalEps = matrix representing the discrete grid of dielectric
    % constants for the simulation
    % featureDims is an array with three values (edge, radius, and triangle
    % a)
    % default is featureDim = celldimx/2 
    epCell = cell(xCells, yCells);
    cellIndices = zeros(xCells, yCells); %%indexes the type of dielectric geometry in the metaatom
    for i = 1:xCells
       for j = 1:yCells
          epsProto = ones(celldimx, celldimy);
          selector = randi([1 3], 1,1);

          cellIndices(i,j) = selector;
          if (selector == 1)
             epsProto = createSquareRod(epsProto,featureDims(1), dielConst); 
          elseif (selector == 2)
             epsProto = createCircularRod(epsProto, featureDims(2), dielConst); 
             
          else
             epsProto = createTriangularRod(epsProto, featureDims(3),dielConst);
          end
          dimCell = size(epsProto);
          %% add the single layer separator on right and bottom 
          epsProto = [epsProto, ones(dimCell(1),1)];
          epsProto = [epsProto; ones(1,dimCell(2)+1)];
          
          epCell{i,j} = epsProto; 
       end
    end
    
    eps_dielect = cell2mat(epCell);
    dim = size(eps_dielect);
    eps_dielect = [ones(dim(1),1), eps_dielect];
    eps_dielect = [ones(1,dim(2)+1); eps_dielect];
    dim = size(eps_dielect);
    
    %% finalEps except we add in the PML dimensions onto it as well
    finalEps = ones(dim(1)+ 2*Npml(1), dim(2) + 2*Npml(2));
    finalEps(Npml(1)+1:end-Npml(1), Npml(2)+1: end-Npml(2)) = eps_dielect;
    
end