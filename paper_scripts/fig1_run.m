hFig = sdd_gui

print('Fig1B.eps','-depsc','-r300');
print('Fig1B.png','-dpng','-r300');

% ax = gca;
exportgraphics(hFig,'Fig1B.pdf','Resolution',300)

print('Fig1B.jpg','-r300');


%%

% test SDD_GUI without any GUI.
sets = Core.Default.read_default_sets('sdd_fig3.txt',0);

sets.ratlim = 0.3;
sets.elim = 0.8;
sets.extractionMethod = 2;
% sets.folder = '/export/scratch/albertas/data_temp/DOTS/031523 h202 data/testFig3/image for slide 7.czi';
% sets.folder = 'C:\Users\Lenovo\test\image for slide 7.czi'
%%
    % create tabbed figure
    hFig = figure('Name', ['SDD-dots GUI v'], ...
        'Units', 'normalized', ...
        'OuterPosition', [0 0 0.3 0.3], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none', ...
        'ToolBar', 'figure' ...
    );

    hPanel = uipanel('Parent', hFig);
    h = uitabgroup('Parent',hPanel);
    t1 = uitab(h, 'title', 'SDD');
    tsHCC = uitabgroup('Parent',t1);
%     hPanelImport = uitab(tsHCC, 'title', 'Dot import tab');

[output,hPanelResult,images,movies,barcodes] = sdd_process_folder(sets.folder , sets,tsHCC);

hPanelResult{1}.Children.SelectedTab  = hPanelResult{1}.Children.Children(2);

% hFig.Position = [100 100 540 400];
% grid on
% axis on
ax = gca;
ax.XLim = [450 660];
ax.YLim = [400 470];
% delete(nexttile(2))
% I2 = imcrop(gcf,[1 100 1 100]);

% gcf(hPanelResult{1}.Children.SelectedTab)

% J = imcrop(hFig,[60 40 100 90]);
set(gcf, 'Color', 'w')

A = gcf;
A.XDisplay;

hold on
nPixels = 1e3/sets.pxnm;
% x = [ax.XLim(2)-nPixels-10 ax.XLim(2)+nPixels-nPixels-10];
x = [ax.XLim(1)+10 ax.XLim(1)+nPixels+10];

y = [ax.YLim(2)-5 ,ax.YLim(2)-5 ];
plot(x,y,'Linewidth',2,'Color','white')
% text(x(1)-2,y(1)-5,'1$\mu$m','FontName','times', 'Fontsize',12,'FontWeight','bold','Color','white','Interpreter','latex')

% text(0,0.05,'10 mM','Fontsize',10,'Color','white','Units','normalized','Interpreter','latex')
set(gcf, 'Color', 'w')

%% Select 14

% ft = 'Times';
% fsz = 10;        
% %%%%%%%%%%%%%%%
% % Your Figure
% %%%%%%%%%%%%%%%
% set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
% set(gca,'FontSize', fsz, 'FontName', ft)



print(A,'Fig1C.png','-dpng','-r300');

%% 14:
sets.pixelSize = sets.pxnm;
Core.AnalysisPlot.detailed_analysis_plot(movies,barcodes,sets,14)
% ax = gca;
% ax.XLim = [450 660];
% ax.YLim = [400 470];

    print('Fig1D.png','-dpng','-r300');

%% Full figure in single run:
%%
f=figure('Position',[1 1 700 800]), % ,'PaperUnits','centimeters'
t = tiledlayout(4,2,'TileSpacing','tight','Padding','tight')
ax0 = nexttile
%     figure(2 + (imageNumber - 1) * 5)
%     axes(tiles.molScores);
    h1 = histogram(log(output{1}.meh(output{1}.meh > 1)));
%     title([imageName, ' edge scores'])

      mess = 'Automated threshold';
      line([log(output{1}.molScoreLim) log(output{1}.molScoreLim)], [0 max(h1.Values)], 'LineStyle', '--', 'Color', 'red', 'LineWidth', 2)

    text(1.1 * log(output{1}.molScoreLim), 2/3 * max(h1.Values), mess, 'FontSize', 12,'FontName','Times')
%     hold off
%     if sets.autoThreshBarsold==0&&sets.autoThreshBarsx`
%         axes(tiles.dotScoresFilt);
%         histogram(allScores)
%         title([imageName, ' edge score background histogram'])
%     end
nexttile
axis off
% daspect(ax0,[1 1 1]);  % <---- move to after-plot
% pbaspect([1 1 1])


ax1 = nexttile([1 2]);
hold on

xlims = [400 742];
ylims = [400 461];
% imagesc(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2)));colormap gray;
imshow(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2)),[min(images{1}.registeredIm{1}(:)) max(images{1}.registeredIm{1}(:))],'InitialMagnification','fit')
% imshow(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2)),[min(images{1}.registeredIm{1}(:)) max(images{1}.registeredIm{1}(:))],'InitialMagnification','fit')

% imshow(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2)),'InitialMagnification','fit')%;colormap gray;

set(gca,'YDir','normal')
axis off
% axis equal
% daspect(ax1,[1 1 1]);  % <---- move to after-plot
% pbaspect([1 0.22 1])

% imshow(mat2gray(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2))), 'InitialMagnification', 'fit');
% nexttile([1 2])
% hold on

% xlims = [450 660];
% ylims = [400 470];
% imshow(mat2gray(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2))), 'InitialMagnification', 'fit');
% nexttile([1 2])
% hold on
% 
% xlims = [450 660];
% ylims = [400 470];
% imshow(mat2gray(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2))), 'InitialMagnification', 'fit');
% 
% figure
% imagesc(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2)));colormap gray;
%     axis equal


ids = [10 14]
for idx=ids
    plot(output{1}.trueedge{idx}(:, 2)-xlims(1)+1, output{1}.trueedge{idx}(:, 1)-ylims(1)+1,'LineWidth',1);
end


hold on
for molIdx = 1:numel(ids)
    %     curIdx = barcodes.idx(molIdx);
    str = sprintf('Mol. %i',molIdx);
    text(-5+min(output{1}.trueedge{ids(molIdx)}(:, 2)-xlims(1)+1),5+max(output{1}.trueedge{ids(molIdx)}(:, 1)-ylims(1)+1),str,'Color','white','FontName','times', 'Fontsize',12, 'Clipping','on');
end


nPixels = 1e3/sets.pxnm;
% x = [ax.XLim(2)-nPixels-10 ax.XLim(2)+nPixels-nPixels-10];
x = [10 nPixels+10];
y = [5 , 5 ];

% y = [diff(ylims)-5 , diff(ylims)-5 ];
plot(x,y,'Linewidth',2,'Color','white')
% text(x(1)-2,y(1)-5,'1$\mu$m','FontName','times', 'Fontsize',12,'FontWeight','bold','Color','white','Interpreter','latex')

% text(0,0.05,'10 mM','Fontsize',10,'Color','white','Units','normalized','Interpreter','latex')

% nexttile
%
sets.pixelSize = sets.pxnm;
Core.AnalysisPlot.detailed_analysis_plot_intensity(movies,barcodes,sets,14)
% ax = gca;
set(gcf, 'Color', 'w')

size(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2)),2)/ size(images{1}.registeredIm{1}(ylims(1):ylims(2),xlims(1):xlims(2)),1)
 size(movies.molM{14},2)/size(movies.molM{14},1)
% 
% nexttile
% imshow(movies.molM{14},[min(movies.molM{14}(:)) max(movies.molM{14}(:))],'InitialMagnification','fit')

% ax = gca;
% copygraphics(ax,'Resolution',300)
    print('Fig1BD.png','-dpng','-r400');
