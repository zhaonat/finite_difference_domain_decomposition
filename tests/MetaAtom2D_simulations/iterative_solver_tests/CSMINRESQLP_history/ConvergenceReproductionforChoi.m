%% Reproduction of Convergence with CSMINRESQLP
close all
clear all

%% load matrices
load('SystemMatrices.mat')

%% perform iterations Analysis
tol = 1e-14;
maxitun = length(b);
maxitred = length(bmod);

%% solve
M = []; %preconditioner
tic
[x1,flag1,iter1,Miter,QLPiter,relres1,relAres,...
      Anorm,Acond1,xnorm,Axnorm,resvec1,Aresvec]  =...
      csminresqlp(A,b,tol, length(b), M);
toc
tic
[x2,flag2,iter2,Miter2,QLPiter2,relres2,relAres2,...
      Anorm2,Acond2,xnorm2,Axnorm2,resvec2,Aresvec2]  = ...
      csminresqlp(Aschur,bmod,tol, 10*length(bmod),M);
toc

%% Visualize Field Solutions
Hz = reshape(x1, length(x1)^.5, length(x1)^.5);
sdf = figure;
subplot(2,1,1)
imagesc(abs(Hz));

f = figure()
semilogy(resvec1/norm(b), 'linewidth', 1.5)
hold on
semilogy(resvec2/norm(bmod), 'linewidth', 1.5)
%title('Convergence History for Reduced and Unreduced Solves')
xlabel('Iteration Number')
ylabel('log_{10}(relative residual)')
legend('unreduced', 'reduced')

 
    



