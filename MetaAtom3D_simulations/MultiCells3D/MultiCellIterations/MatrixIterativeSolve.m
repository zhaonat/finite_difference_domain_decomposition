%% iterative matrix solve

function soln = MatrixIterativeSolve(A,B)
    d = size(B);
    numCols = d(2);
    soln = [];
    for i =1:numCols
        tic
        y = qmr(A,B(:,i),1e-8, d(1));
        toc
        soln = [soln,y];
    end

end