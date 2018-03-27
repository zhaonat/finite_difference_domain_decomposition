%% coordinate indexing multicell analysis
Nx = 500; Ny = 500;
xCells = 10; yCells = 10;
SCS = Nx/xCells;
tic
[Nindexing, bC, iC] = CoordinateIndexing2DMultiCell(Nx, Ny, xCells, yCells);
toc

for i = 0:xCells*yCells-1
   plot(Nindexing(i*SCS^2+1:(i+1)*SCS^2,1), Nindexing(i*SCS^2+1:(i+1)*SCS^2,2),...
       '.', 'markersize', 15, 'color', rand(3,1))
   hold on;
   
end