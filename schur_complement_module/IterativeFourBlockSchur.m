

%% the operation inside here which calculates the schur complement
%% is parallelizable!

function [Aschur, bmod, App, Avvcell, Apvcell, Avpcell, bvCell] = ...
    IterativeFourBlockSchur(SymA,SymB, hpart,vpart, xCells, yCells, N)
 Nx = N(1); Ny = N(2);
CDx = Nx/xCells; CDy = Ny/yCells;
totCells = xCells*yCells;


%uppermost block = boundary cell partition
    % lower blocks = interior cell partition
    Abound = SymA(1:hpart, 1:vpart);
    

    Aint = SymA(hpart+1:end, vpart+1:end);
    Aib = SymA(hpart+1:end,1:vpart);
    Abi = SymA(1:hpart, vpart+1:end);
    
    bBound = SymB(1:hpart);
    bint = SymB(hpart+1:end);
   
    App = Abound;
    Avvcell = mat2cell(Aint, (CDx-2)^2*ones(totCells,1), (CDy-2)^2*ones(totCells,1));
    Apvcell = mat2cell(Abi, hpart, (CDx-2)^2*ones(totCells,1));
    Avpcell = mat2cell(Aib, (CDy-2)^2*ones(totCells,1), vpart);
    bvCell = mat2cell(bint, (CDx-2)^2*ones(totCells,1), 1);
    bp = bBound;
    Aschur = App;
    bmod = bp;
    to = cputime;
    
    parfor i = 1:xCells*yCells
        
        A2 = Avvcell{i,i};
        Apv = Apvcell{1, i};
        Avp = Avpcell{i, 1};
        Avv = A2;
        bp = bBound;
        bv =bvCell{i,1};
        Comp = Apv*(Avv\Avp);
        Aschur = Aschur - Comp; 
        bmod = bmod - Apv*(Avv\bv);
       
    end
    
end    