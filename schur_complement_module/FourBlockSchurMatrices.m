%% 2x2 block schur complement

function [App, Apv, Avp, Avv, bp, bv, bmod] = ...
    FourBlockSchurMatrices(A,b, hpart, vpart)
    %A = [ 1 2 3 4; 5 6 7 8; 9 1 2 3; 6 5 4 2];
    %b = [ 1 2 3 4];
   
    dimension = size(A);
    App = A(1:hpart, 1:vpart);
    Apv = A(1:hpart, vpart+1:dimension(2));
    Avp = A(hpart+1:dimension(1), 1:vpart);
    Avv = A(hpart+1:dimension(1),vpart+1:dimension(2));
    
    bp = (b(1:hpart));
    bv =(b(hpart+1:dimension(2)));
    bmod = bp -Apv*(Avv\bv);
    % = circshift(Avv, 110) circshifting does not affect condition number
    
    %%reconstruct original solution
%     sol1 = ASchur\bmod;
%     sol2 = (Avv\(bv - Avp*sol1));

end    