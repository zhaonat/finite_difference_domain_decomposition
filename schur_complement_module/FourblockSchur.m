%% 2x2 block schur complement

function [Aschur, bmod, App, Apv, Avp, Avv, bp, bv] = ...
    FourblockSchur(A,b, hpart)
    %A = [ 1 2 3 4; 5 6 7 8; 9 1 2 3; 6 5 4 2];
    %b = [ 1 2 3 4];
    vpart = hpart; 
    dimension = size(A);
    App = A(1:hpart, 1:vpart);
    Apv = A(1:hpart, vpart+1:dimension(2));
    Avp = A(hpart+1:dimension(1), 1:vpart);
    Avv = A(hpart+1:dimension(1),vpart+1:dimension(2));
    disp(strcat('size of AVV=',num2str(length(Avv))));
    %[L,U] = lu(Avv);
    bp = (b(1:hpart));
    bv =(b(hpart+1:dimension(2)));
    Aschur = App - Apv*(Avv\Avp);
    bmod = bp-Apv*(Avv\bv);
%     [L,D,P] = ldl(Avv);
%     ASchur = App - Apv*(P.'\(L.'\(D\(L\(P\Avp))))); %% this is a potentially costly operation
%     %%as you have the inverse of Avv, a large square matrix
%     bmod = bp -Apv*(P.'\(L.'\(D\(L\(P\bv)))));
    
    %%reconstruct original solution
%     sol1 = ASchur\bmod;
%     sol2 = (Avv\(bv - Avp*sol1));

end    