%% Dielectric MultiCell Formatter
%% generates a layered grid of dielectrics (need it for consistency in paper)
%% in the domain and the protoypical cell

function [finalEps, cellIndices] = ...
    UniformLayerCellDielectricSingleLayerSep(xCells, yCells,celldimx, celldimy, Npml, dielConst)

    %% Parameter Descriptions
    % xCells - # of metaatoms in x direction
    % yCells - # of metaatoms in y direction
    % cellindices - a dictionary of what dielectric inclusions is
    %     contained in which meatatom; this is a matrix with the dimension of
    %     the metatom grid: xcells*ycells;
    % finalEps = matrix representing the discrete grid of dielectric
    % constants for the simulation
    
    epCell = cell(xCells, yCells);
    cellIndices = zeros(xCells, yCells); %%indexes the type of dielectric geometry in the metaatom
    featureDim = celldimx/2;
    for i = 1:xCells
       for j = 1:yCells
          epsProto = ones(celldimx, celldimy);
          
          %% selector index dictates what cell gets what geometry,
           % this is essentially the god of the dielectric generator
          selector = 0;
          if(i < ceil(xCells/3) == 1)
             selector = 2;
          elseif(i >= ceil(xCells/3) && i <ceil(2*(xCells/3)) )
              selector = 3;
          else
              selector = 1;
          end

          cellIndices(i,j) = selector;
          if (selector == 1)
             epsProto = createSquareRod(epsProto,featureDim, dielConst); 
          elseif (selector == 2)
             epsProto = createCircularRod(epsProto, featureDim/2, dielConst); 
             
          else
             epsProto = createTriangularRod(epsProto, featureDim/3,dielConst);
          end
          dim = size(epsProto);
          epsProto = [epsProto, ones(dim(1),1)];
          epsProto = [epsProto; ones(1,dim(2)+1)];
          
          epCell{i,j} = epsProto; 
       end
    end
    
    eps_dielect = cell2mat(epCell);
    dim = size(eps_dielect);
    eps_dielect = [ones(dim(1),1), eps_dielect];
    eps_dielect = [ones(1,dim(1)+1); eps_dielect];
    dim = size(eps_dielect);
    
    %% finalEps except we add in the PML dimensions onto it as well
    finalEps = ones(dim(1)+ 2*Npml(1), dim(2) + 2*Npml(2));
    finalEps(Npml(1)+1:end-Npml(1), Npml(2)+1: end-Npml(2)) = eps_dielect;
    
end