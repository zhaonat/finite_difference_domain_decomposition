function k = multiCellVisualization(recSolIntArray, xrange, yrange)
    k = figure;
    dim = size(recSolIntArray);
    count = 1;
    for i = 1:length(recSolIntArray)
        for j = 1:length(recSolIntArray)

            subplot(dim(1), dim(2), count)
            imagesc(abs(recSolIntArray{j,i}))
            %visabs(recSolIntArray{i,j}, xrange/dim(1)-1, yrange/dim(2)-1);
            count = count+1;
        end
    end
end