
function J = PlaneWave(N, k, omega, c)

    kx = k(1);
    ky = k(2);
    xa = 1:N(1);
    ya = 1:N(2);
    [Y,X] = meshgrid(ya,xa);
    k0 = omega/c;
    J = exp(-1i*(k0*(kx*X+ky*Y)));
    
end