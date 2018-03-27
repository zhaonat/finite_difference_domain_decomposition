%% Dielectric MultiCell Formatter
%% a function which generates a 2D eps_dielectric matrix given the number of cells
%% in the domain and the protoypical cell

function [finalEps, cellIndices] = multiRandomCellDielectric(xCells, yCells,...
    celldimx, celldimy, Npml, epsilon, featureDims)
    epCell = cell(xCells, yCells);
    cellIndices = zeros(xCells, yCells);
    %featureDim = celldimx/2;
    for i = 1:xCells
       for j = 1:yCells
          epsProto = ones(celldimx, celldimy);
          selector = randi([1 3], 1,1);
          cellIndices(i,j) = selector;
          if (selector == 1)
             epsProto = createSquareRod(epsProto,featureDims(1), epsilon); 
          elseif (selector == 2)
             epsProto = createCircularRod(epsProto, featureDims(2), epsilon); 
             
          else
             epsProto = createTriangularRod(epsProto, featureDims(3),epsilon);
          end
          
          epCell{i,j} = epsProto; 
       end
    end
    
    
    eps_dielect = cell2mat(epCell);
    dim = size(eps_dielect);

    finalEps = ones(dim(1)+ 2*Npml(1), dim(2) + 2*Npml(2));

    finalEps(Npml(1)+1:end-Npml(1), Npml(2)+1: end-Npml(2)) = eps_dielect;
end