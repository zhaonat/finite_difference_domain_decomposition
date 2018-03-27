%% dielectric multicell analysis

xCells = 2; yCells = 2;
epsPrototype = ones(5,5)
Npml = [2 2];
celldimx = 60; celldimy = 60;
dielConst = 12;
featureDims = [40, 20, 20];
[finalEps, cellIndices] = ...
    multiRandomCellDielectricSingleLayerSep(xCells, yCells,celldimx, ...
    celldimy, Npml, dielConst, featureDims)

figure;
imagesc(finalEps);