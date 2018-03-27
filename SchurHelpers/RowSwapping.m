%% Matrix Swapping Function
function swapped = RowSwapping(A,b)

AUG = [A b];
interiorPoints = [];
boundaryPoints = [];
M= N(1)*N(2);

for i = 1:N(1)*N(2) %x index
    if(i == 1 || i <= N(1) || mod(i, N(1)) == 0 || (i> M-N(1) && i<M))
       boundaryPoints = sparse(vertcat(boundaryPoints, AUG(i,:))); 
    else
       interiorPoints = sparse(vertcat(interiorPoints, AUG(i,:)));
    end
   
end
swapped = [interiorPoints; boundaryPoints];


end