    N = [8 8 8];
     indexorder = transpose(1:N(1)*N(2)*N(3));
    
    %% Call IndexPermutation Function to achieve reordering
    permutedIndices = IndexPermutation3D_3Field(N);
    permutedIndicesSingle = IndexPermutation3D(N);
    %% Execute a Symmetry Preserving Column Row Permutation Combination
    Q = speye(3*N(1)*N(2)*N(3));
    for i = 1:3*N(1)*N(2)*N(3)
       Q(i,i) = 0;
       indexshift = permutedIndices(i);
       Q(i,indexshift) = 1;
    end
    Q2 = speye(N(1)*N(2)*N(3));
    for i = 1:N(1)*N(2)*N(3)
       Q2(i,i) = 0;
       indexshift = permutedIndicesSingle(i);
       Q2(i,indexshift) = 1;
    end
    Q3 = blkdiag(Q2,Q2,Q2);
    %%test permutation matrix should permute index order to permuted indices
    indexorder = [indexorder;indexorder;indexorder];
    y = Q*indexorder;
    
    %% Transform Equations
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    symB = Q*b;
    
    SymA2=Q3*A*transpose(Q3);
    SymB2 = Q3*b;