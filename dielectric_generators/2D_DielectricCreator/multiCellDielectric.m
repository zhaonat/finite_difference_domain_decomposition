%% Dielectric MultiCell Formatter
%% a function which generates a 2D eps_dielectric matrix given the number of cells
%% in the domain and the protoypical cell

function finalEps = multiCellDielectric(xCells, yCells, epsPrototype, Npml)
    epCell = cell(xCells, yCells);
    for i = 1:xCells
       for j = 1:yCells
          epCell{i,j} = epsPrototype; 
       end
    end
    eps_dielect = cell2mat(epCell);
    dim = size(eps_dielect);

    finalEps = ones(dim(1)+ 2*Npml(1), dim(2) + 2*Npml(2));

    finalEps(Npml(1)+1:end-Npml(1), Npml(2)+1: end-Npml(2)) = eps_dielect;
end