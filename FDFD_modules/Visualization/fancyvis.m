function [h, Outline] = fancyvis(Hz, eps_air, xCells, SingleCellSize, N, Npml,...
    outline,fontsize, colorbarStat)
    
    Nx = N(1); Ny = N(2);
    % construct structure outline
    Outline = rot90(fliplr(edge(abs(eps_air),'Sobel')));
    index = [];
    for i = 1:Nx
       for j = 1:Ny
          if(Outline(i,j) == 1)
             index = [index;i,j]; 
          end
       end
    end

    h = imagesc(real(Hz));
    pbaspect([1 1 1])
    cmax = max(abs(Hz(:)));
    caxis([-cmax, cmax]);
    colormap('b2r')

    pbaspect([1 1 1])
    set(gca,'Ydir','Normal')
    set(gca, 'FontSize', fontsize)

    c = 0*[1,1,1];
    xlabel('\mu m', 'Color',c)
    ylabel('\mu m', 'Color', c)
    ax.XColor = 0.0*[1,1,1];
    ax.YColor = 0.0*[1,1,1];
    
    tickGrid = [];
    if(sum(Npml) >0)
        for i = 0:2:xCells
           tickGrid = [tickGrid, Npml(1)+i*SingleCellSize+i+1]; 
        end
    else
        for i = 0:2:xCells
           tickGrid = [tickGrid, i*SingleCellSize+i+1]; 
        end
    end
    if(sum(Npml) > 0)
        tickGrid = [tickGrid, tickGrid(end)+Npml(2)];
    end
    tickGrid
    set(gca,'Xtick',tickGrid)
    set(gca,'Ytick',tickGrid)
    xticklabels([0:4:2*xCells]);
    yticklabels([0:4:2*xCells]);
    hold on;
    if(nnz(Outline) > 0)
        plot(index(:,1), index(:,2), '.k', 'markersize', outline);
    end
    %% ADD PML Boundary
    line([Npml(1) Npml(2)], [0, N(1)], 'linewidth', 1.25, 'Color','black','LineStyle','--')
    line([0, N(1)], [Npml(1) Npml(2)],  'linewidth', 1.25, 'Color','black','LineStyle','--')
    line([N(1)-Npml(1) N(1)-Npml(2)], [0, N(1)], 'linewidth', 1.25, 'Color','black','LineStyle','--')
    line([0, N(1)], [N(1)-Npml(1) N(1)-Npml(2)], 'linewidth', 1.25, 'Color','black','LineStyle','--')
    set(gca, 'LineWidth',1.25)
    cmax = max(max(abs(Hz)));
    cmax = full(cmax);
    
    
    %% colorbar
    if(colorbarStat == 1)
        hcb = colorbar;
        set(hcb, 'YTick', [-cmax,0, cmax])
        set(hcb,'TickLabels',{'-max', '0', 'max'})
    end
end