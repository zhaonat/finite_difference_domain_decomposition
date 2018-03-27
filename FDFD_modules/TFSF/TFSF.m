function [E_TFSF] = TFSF(Q,A,Jz,L0,wvlen,xrange, yrange, Npml)

    %% solve vacuum problem
    N = size(Jz); %JZ is still a MATRIX
    eps_r = ones(N);
    [ A0, ~,b] =...
        solveTM_Matrices(L0,wvlen, xrange, yrange, eps_r, Jz, Npml);
    fsrc = A0\b; %vacuum field distribution of the point source
    
    %% solve TFSF problem

    bprime = (Q*A - A*Q)*fsrc;

    %% solve TFSF problem
    x = A\bprime; %% the tfsf problem isolates ONLY  the scattered fields...?

    E_TFSF = reshape(x, N(1),N(2));
    
end