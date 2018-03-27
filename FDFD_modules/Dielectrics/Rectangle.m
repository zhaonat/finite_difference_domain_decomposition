
function Grid = Rectangle(Grid, r1, r2, epdiel)
    
    Coords = [r1;r2];
    xmin = min(Coords(:,1));
    xmax = max(Coords(:,1));
    ymin = min(Coords(:,2));
    ymax = max(Coords(:,2));
    Grid(ymin:ymax,xmin:xmax) = epdiel;
    
end