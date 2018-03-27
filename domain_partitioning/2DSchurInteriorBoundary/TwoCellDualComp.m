function [SymA SymB, B1Coordinates, D1Coordinates, B2Coordinates, D2Coordinates] ...
    = TwoCellDualComp(A,b, N)
    NaturalOrderingCoordinates = CoordinateIndexing(N(1), N(2));
    B1Coordinates = []; B2Coordinates = []; D1Coordinates = []; D2Coordinates = [];
    permutedIndices = []; B1Indices = []; B2Indices = []; D1Indices = []; D2Indices = [];
    for i = 1:length(NaturalOrderingCoordinates)
       x = NaturalOrderingCoordinates(i,1);
       y = NaturalOrderingCoordinates(i,2);

       if(x <= N(1)/2)
           if(x == 1 || y == 1 || y == N(2) || x == N(1)/2)
               B1Coordinates = [B1Coordinates; x y];
               B1Indices = [B1Indices; i];
           else
               D1Coordinates = [D1Coordinates; x y];
               D1Indices = [D1Indices; i];
           end
       else
           if(x == N(1) || y == 1 || y == N(2) || x == N(1)/2+1)
               B2Coordinates = [B2Coordinates; x y];
               B2Indices = [B2Indices; i];
           else
               D2Coordinates = [D2Coordinates; x y];
               D2Indices = [D2Indices; i];
           end

       end
   
    end
    permutedIndices = [B1Indices; D1Indices; B2Indices; D2Indices]


    %% Execute a Symmetry Preserving Column Row Permutation Combination
    Q = speye(N(1)*N(2));
    for i = 1:N(1)*N(2)
       Q(i,i) = 0;
       indexshift = permutedIndices(i);
       Q(i,indexshift) = 1;
    end

    %%test permutation matrix should permute index order to permuted indices

    %% Transform Equations
    SymA = Q*A*transpose(Q);
    %issymmetric(SymA)
    SymB = Q*b;



end