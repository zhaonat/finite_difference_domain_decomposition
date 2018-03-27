%% input is r an (d,1) vector of the coordinate in the cell
% and N a (d,1) vector of the rectangular domain size

function k = CoordtoNatIndex(r,N)
    if(length(N) == 3) %% 3D case
        k = (r(1)) +(r(2)-1)*N(2)+(r(3)-1)*(N(1)*N(2)) ;
    else %% 2D case
        k = (r(1)) + (r(2)-1)*N(2);
    end
end