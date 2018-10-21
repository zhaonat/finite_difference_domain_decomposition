%% Function to take the multicell boundary solution and reconstruct full solution
%% using the reverse permutation matrix
%% Function INputs
%% Q = permutation matrix
%% xschur = schur solution
%% Avp, bv are the proper subblocks of the whole system

function recField = reconstruct_four_block(xSchur, Avv, Avp, bv, Q)

    yVec = bv - Avp*xSchur;
    %[recSolInt, flagI, relresI, iterI, resvecI] = qmr(Avv,yVec, 1e-8, 1e8);
    recSolInt = Avv\yVec;
    tolSol = [xSchur; recSolInt];
    
    recField = Q\tolSol;
end