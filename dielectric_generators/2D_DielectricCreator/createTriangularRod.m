function eps_proto = createTriangularRod(Grid, a, epdiel)
    %% a is a parameter which determines the 'sidelength'
    dim = size(Grid);
    Nx = dim(1);
    Ny = dim(2);
    midptx = ceil(Nx/2);
    midpty = ceil(Ny/2);
    h = sqrt(3)*a/2;
    f = sqrt(3)/2;
    for i = 1:Nx
       for j = 1:Ny
           x = i;
           y = j;
           if(y >= midpty - h && x >= midptx - h && y<= -(x)*f +2 *midptx)
               Grid(x,y) = epdiel;
           end
       end
        
    end
    eps_proto = Grid;
end
