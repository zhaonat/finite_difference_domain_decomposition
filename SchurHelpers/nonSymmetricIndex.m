
%% given a matrix, returns indices where the matrix is NOT symmetric

%% given A
function [wrongIndices] = nonSymmetricIndex(A)
    [xcoor, ycoor, val] = find(A);
    wrongIndices = []; 
    for i = 1:length(xcoor)
        x = xcoor(i);
        y = ycoor(i);
        value = val(i);

        xtrans = y;
        ytrans = x;
        for j = i:length(xcoor)
           if(xcoor(j) == xtrans && ycoor(j) == ytrans && xcoor(j) ~= ycoor(j))
%               disp(strcat('x', int2str(xcoor(i))));
%               disp(strcat('y', int2str(ycoor(i))));
              if(abs(val(j)-value) > eps(max(max(A)))*10)
                  wrongIndices = [wrongIndices; xcoor(j) ycoor(j) val(j) value];
              end
           end
        end
    % end of outer loop
    end
end