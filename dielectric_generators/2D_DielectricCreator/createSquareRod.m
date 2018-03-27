
function eps_proto = createSquareRod(Grid, edgeLength, epdiel)
    
    dim = size(Grid);
    Nx = dim(1);
    Ny = dim(2);
    midptx = ceil(Nx/2);
    midpty = ceil(Ny/2);
    Grid(midptx-ceil(edgeLength/2)+1:midptx+ceil(edgeLength/2), ...
        midpty-ceil(edgeLength/2)+1:midpty+ceil(edgeLength/2)) = epdiel;
    eps_proto = Grid;

end