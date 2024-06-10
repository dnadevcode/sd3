function [] = detailed_analysis_plot_intensity(movies,barcodes,sets,idx)

    nmPx = 1/sets.pixelSize*1000;

%     figure,
%     idx = 62;
    mol = movies.molM{idx};
    mol2 = movies.dotM{idx};
    
    % if there are some delid, the correct idx has to be checked
    % (barcodes.idx)

    xF = barcodes.xy{idx}{1};
    yF = barcodes.xy{idx}{2};

    nexttile
    plot(barcodes.expBars{idx}.rawBarcode(~isnan(barcodes.expBars{idx}.rawBarcode)))
    title('YOYO-1 intensity','fontname','Times','Interpreter','latex')
%     xgr = [0:nmPx:length(barcodes.dotBars{idx})];
%     xticks(xgr)
%     xticklabels([arrayfun(@(x) num2str(x),0:1:length(xgr),'un',false)] )
%     xtickangle(0)
%     xlabel('Position  $\mu m$','fontname','Times','Interpreter','latex')
    ylabel('Intensity','fontname','Times','Interpreter','latex')
%     axis off
%     axis equal
%     pbaspect([1 1 1])
h = gca;
    set(gca,'xtick',[])

% h.XAxis.Visible = 'off';
    nexttile
    plot(barcodes.dotBars{idx}(~isnan(barcodes.dotBars{idx})))
    title('Dot intensity','Interpreter','latex')
%     h = gca;
    set(gca,'xtick',[])
%     h.XAxis.Visible = 'off';
%     xgr = [0:nmPx:length(barcodes.dotBars{idx})];
%     xticks(xgr)
%     xticklabels([arrayfun(@(x) num2str(x),0:1:length(xgr),'un',false)] )
%     xtickangle(0)
%     xlabel('Position $\mu m$','fontname','Times','Interpreter','latex')
%     ylabel('Intensity','fontname','Times','Interpreter','latex')
%     nexttile
%     %x = ["Conveg to hull" "Eccentricity" "Length"];
%     bar([movies.stats{idx}.FilledArea/movies.stats{idx}.ConvexArea movies.stats{idx}.Eccentricity length(barcodes.expBars{idx}.rawBarcode)/100])

%     axis off
%     axis equal
%     pbaspect([1 1 1])

    ax2 =nexttile([1 2]);
    imshow(mol2,[min(mol2(:)) max(mol2(:))],'InitialMagnification',100)

    hold on
%     imagesc(mol2)
    set(gca,'YDir','normal')

    hold on
    %             
    plot(yF,xF,'->','MarkerSize',3, 'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
    % %             
    colormap(gray)
    pos = barcodes.dots{idx}.locations+barcodes.dots{idx}.leftOffset;
%     pos
    plot(barcodes.xy{idx}{2}(pos),barcodes.xy{idx}{1}(pos),'greeno','MarkerSize',20)
%     daspect(axA,[1 1 1]);  % <---- move to after-plot
%     pbaspect(axA,[1 1 1]); % <---- move to after-plot
    ylim([1 size(mol,1)])
    xlim([1 size(mol,2)])
    xgr = [0:nmPx:length(barcodes.dotBars{idx})];
    xticks(xgr)
    xticklabels([arrayfun(@(x) num2str(x),0:1:length(xgr),'un',false)] )
    xtickangle(0)
%     xlabel('Position (micron)','fontname','Times')
%     axis equal

    x = [5 nmPx+5];
%     y = [size(mol,1)-2 size(mol,1)-2];
    y = [2 2];

    plot(x,y,'Linewidth',2,'Color','white')

    axis off

% axis equal
% daspect(ax2,[1 1 1]);  % <---- move to after-plot
% pbaspect([4 0.4 1])
%     ft = 'Times';
%     fsz = 10;        
   %%%%%%%%%%%%%%%
    % Your Figure
    %%%%%%%%%%%%%%%
%     set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
%     set(gca,'FontSize', fsz, 'FontName', ft)


end