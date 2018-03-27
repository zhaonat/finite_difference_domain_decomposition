%C is a 2x2 cell of matrix subblocks
%function does a schur complement of lower right to upper left
% b is a 2x1 cell containing two vectors
function [Aschur,bmod,L,U] = FourblockSchurCell(C, b)
    %do ilu factorization of C{2,2}
    setup.type = 'crout';
    setup.milu = 'row';
    setup.droptol = 0.0;
    [L,U,P] = ilu(sparse(C{2,2}), setup);
    C{2,2} = L*U;
    Aschur = C{1,1} - C{1,2}*(U\(L\C{2,1}));
    a = b{1}; b= b{2};
    bmod = a - C{1,2}*(U\(L\b));
end