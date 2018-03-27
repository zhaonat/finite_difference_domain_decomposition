function [x,y] = NatIndexToCoord(k,Nx, Ny)
    y = mod(k,Ny);
    if(y==0)
       y = Ny; 
    end
    x = ceil(k/Nx);
end