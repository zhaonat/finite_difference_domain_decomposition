%% iterative matrix solve

function soln = MatrixIterativeSolve(A,B)
    d = size(B);
    numCols = d(2);
    soln = [];
    tic
    for i =1:numCols
        
        [y,~,~] = qmr(A,B(:,i),1e-8, d(1));
        
        soln = [soln,y];
    end
    toc
end