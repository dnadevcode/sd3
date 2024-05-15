function [] = detailed_analysis_plot(movies,barcodes,sets,idx)

% 
% indexes = {62};
%     for ii=1:length(indexes)
%           idx=find(barcodes.idx==indexes{ii})
%           figure,imagesc(movies.dotM{barcodes.idx(idx)})
%           hold on
%           plot(barcodes.xy{barcodes.idx(idx)}{2},barcodes.xy{barcodes.idx(idx)}{1},'redx')
%           pos = barcodes.dots{idx}.locations+barcodes.nanid(idx);
%           plot(barcodes.xy{idx}{2}(pos),barcodes.xy{idx}{1}(pos),'blackx')
%     %             plot(-barcodes.lineParams{idx}(1)*(1:size(movies.molM{idx},2))+barcodes.lineParams{idx}(2),'redx')
% 
% %           plot(-barcodes.lineParams{idx}(1)*(1:size(movies.molM{idx},2))+barcodes.lineParams{idx}(2),'redx')
%     % %     %   
%         colormap(gray)
% 
% 
% end

    nmPx = 1/sets.pixelSize*1000;

    figure,
%     idx = 62;
    mol = movies.molM{idx};
    mol2 = movies.dotM{idx};
    
    % if there are some delid, the correct idx has to be checked
    % (barcodes.idx)

    xF = barcodes.xy{idx}{1};
    yF = barcodes.xy{idx}{2};
    
    tiledlayout(2,3);axA =nexttile;
    imagesc(mol)
    hold on
    axis equal
    ylim([1 size(mol,1)])
%     daspect(axA,[1 1 1]);  % <---- move to after-plot
%     pbaspect(axA,[1 1 1]); % <---- move to after-plot

    axA =nexttile;

    imagesc(mol2)
    hold on
    %             
    plot(yF,xF,'->','MarkerSize',3, 'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
    % %             
    colormap(gray)
    pos = barcodes.dots{idx}.locations+barcodes.dots{idx}.leftOffset;
%     pos
    plot(barcodes.xy{idx}{2}(pos),barcodes.xy{idx}{1}(pos),'greeno','MarkerSize',20)
    axis equal
%     daspect(axA,[1 1 1]);  % <---- move to after-plot
%     pbaspect(axA,[1 1 1]); % <---- move to after-plot
    ylim([1 size(mol,1)])
    xgr = [0:nmPx:length(barcodes.dotBars{idx})];
    xticks(xgr)
    xticklabels([arrayfun(@(x) num2str(x),0:1:length(xgr),'un',false)] )
    xtickangle(0)
    xlabel('Position (micron)','fontname','Times')
    axA = nexttile;
    imagesc(movies.bwM{idx})
    axis equal
%     daspect(axA,[1 1 1]);  % <---- move to after-plot
%     pbaspect(axA,[1 1 1]); % <---- move to after-plot
    ylim([1 size(mol,1)])

    nexttile
    plot(barcodes.expBars{idx}.rawBarcode(~isnan(barcodes.expBars{idx}.rawBarcode)))
    title('YOYO-1 intensity plot','fontname','Times')
    xgr = [0:nmPx:length(barcodes.dotBars{idx})];
    xticks(xgr)
    xticklabels([arrayfun(@(x) num2str(x),0:1:length(xgr),'un',false)] )
    xtickangle(0)
    xlabel('Position (micron)','fontname','Times')
    ylabel('Intensity','fontname','Times')

    nexttile
    plot(barcodes.dotBars{idx}(~isnan(barcodes.dotBars{idx})))
    title('Dot intensity plot')

    xgr = [0:nmPx:length(barcodes.dotBars{idx})];
    xticks(xgr)
    xticklabels([arrayfun(@(x) num2str(x),0:1:length(xgr),'un',false)] )
    xtickangle(0)
    xlabel('Position (micron)','fontname','Times')
    ylabel('Intensity','fontname','Times')
    nexttile
    %x = ["Conveg to hull" "Eccentricity" "Length"];
    bar([movies.stats{idx}.FilledArea/movies.stats{idx}.ConvexArea movies.stats{idx}.Eccentricity length(barcodes.expBars{idx}.rawBarcode)/100])


    ft = 'Times';
    fsz = 10;        
    %%%%%%%%%%%%%%%
    % Your Figure
    %%%%%%%%%%%%%%%
    set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
    set(gca,'FontSize', fsz, 'FontName', ft)


end