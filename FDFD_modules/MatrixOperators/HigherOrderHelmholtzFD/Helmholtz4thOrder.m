%% source is from I.Singer, E. Turkel Tel Aviv University

function A = Helmholtz4thOrder(omega,N, dL)
    h = dL(1); %assume equal spacing in y and x
    k = omega;
    Laplacian = Laplacian9pointNN(dL,N);
    
    %% modify the diagonal coefficient
    M = prod(N);
    diagCorrection = spdiags((k*h)^2*(1-(k*h)^2/12)*ones(M,1),0,M,M);
    diagCorrection = diagCorrection/h^2;
    
    HomogeneousHelmholtz = Laplacian +diagCorrection;
    A = HomogeneousHelmholtz;

end