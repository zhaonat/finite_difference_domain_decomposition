%% numerical dispersion analysis for FDFD
% taken from taflove
close all

dL = [0.1, 0.1];
dx = dL(1); dy = dL(2);
delta = dx;
figure; 
c = 1;
Ns = [];
for delta = 2*[0.05, 0.1]
    lambda_0 = 1;
    dt = 0.01; %% dt displaces the curves up and down

    %% newton iteration for k icount

    phi = linspace(0,pi/2, 100);

    S = c*dt/delta;
    N_lambda = lambda_0/delta;
    A = delta*cos(phi)/2;
    B = delta*sin(phi)/2;
    C = (1/S^2)*sin(pi*S/N_lambda)^2;

    k_count = 2*pi*ones(length(phi),1).';

    for i = 1:5
       update = (sin(A.*k_count).^2 + sin(B.*k_count).^2 - C)./...
           (A.*sin(2*A.*k_count) + B.*sin(2*B.*k_count));
       k_count = k_count - update; 
    end


    plot(phi, 2*pi./k_count);
    hold on;
    
    %% approximate reflection coefficient
    Ns = [Ns; k_count];
    
end

%assume  k_count is proportional to the index of refraction
R = (Ns(2,:) - Ns(1,:))./(Ns(2,:)+Ns(1,:));
figure; 
plot(R);
