
function eps_proto = createCircularRod(Grid, rad, epdiel)
    
    dim = size(Grid);
    Nx = dim(1);
    Ny = dim(2);
    midptx = ceil(Nx/2);
    midpty = ceil(Ny/2);
    for i = 1:Nx
       for j = 1:Ny
           if((i-midptx)^2 + (j-midpty)^2 <= rad^2)
              Grid(i,j)= epdiel; 
           end
       end
        
    end
    
    eps_proto = Grid;
end
