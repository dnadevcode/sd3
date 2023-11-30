function [] = fig_manual_vs_sdd()

% manual counting vs sdd




xlsFiles = dir(['/export/scratch/albertas/data_temp/DOTS/data_for_figures/Manual/','*.xlsx']);
for i=1:length(xlsFiles)
    dataSSDmanual{i} = readcell(fullfile(xlsFiles(i).folder,xlsFiles(i).name));
    dotsSDD = cell2mat(dataSSDmanual{i}(3:end,1));
    dotsManual = cell2mat(dataSSDmanual{i}(3:end,2));
end

% 
% figure,plot(dotsSDD,dotsManual,'x')
% 
x = dotsSDD';
y =  dotsManual';
p = ranksum(x,y)
% %Engine
[uxy, jnk, idx] = unique([x.',y.'],'rows');
szscale = histcounts(idx,min(unique(idx))-0.5:max(unique(idx))+0.5);
%Plot Scale of 25 and stars
f=figure
tiledlayout(1,2)
nexttile
hold on
ax = scatter(uxy(:,1),uxy(:,2),'c','blue','filled','sizedata',szscale*3)
ax.MarkerEdgeColor = [0.85  0.32 0.098];
ax.MarkerFaceColor = [0.85  0.32 0.098];
xlabel('SDD Number of Dots')
ylabel('Manual Number of Dots')

mdl = fitlm(x,y)
% text(1,max(y)-1,['f(x) = ',num2str(mdl.Coefficients.Estimate(2)),'x','+',num2str(mdl.Coefficients.Estimate(1))])
% text(1,max(y)-2,['R^2 =',num2str(mdl.Rsquared.Ordinary)])

plot(min(x):max(x),mdl.Coefficients.Estimate(2)*(min(x):max(x))+mdl.Coefficients.Estimate(1))
lgd=legend({['$R^2$ =',num2str(mdl.Rsquared.Ordinary)],['f(x) = ',num2str(mdl.Coefficients.Estimate(2)),'x','+',num2str(mdl.Coefficients.Estimate(1))]},'Interpreter','latex')
lgd.Location ='southoutside';
% 
% figure
% tiledlayout(2,2)
% nexttile
% imagesc(A1)
% nexttile
% imagesc(B);colormap(gray)
% title('Mol 145')
% nexttile
% scatter(lengthsM(keepRows),calculatedLengths(keepRows))
% xlabel('Length SDD')
% ylabel('Length tiff16')
% nexttile
% ax=scatter(uxy(:,1),uxy(:,2),'filled','c','sizedata',szscale*25)
% ax.MarkerEdgeColor = [0.85  0.32 0.098];
% ax.MarkerFaceColor = [0.85  0.32 0.098];
% xlabel('Num dots SDD')
% ylabel('Num dots tiff16')

% nexttile

% [output,hPanelResult,images,movies,barcodes]  = sdd_script('sdd_fig_obed.txt');

    
% f=figure
% tiledlayout(2,3,'TileSpacing','compact');
% nexttile([2 2])
% imagesc(images{1}.dotIm);colormap(gray)
% hold on
% trueedge = movies.trueedge;
% %     axes(tiles.dotDet);
% hold on
% 
% for k = 1:length(trueedge)
% %       figure(dotFigNum)
%   plot(trueedge{k}(:, 2), trueedge{k}(:, 1));
% end
% 
% for molIdx = 1:length(barcodes.dots)
%     curIdx = barcodes.idx(molIdx);
%     str = sprintf('Mol. %i',curIdx);
%     angle = atan(barcodes.lineParams{curIdx}(1));
%     % 
%     vOff = barcodes.boundaries{curIdx}(1);
%     hOff = barcodes.boundaries{curIdx}(3);
%     voffList(molIdx) = vOff;
% end
% 
% hold on
%  for molIdx = 1:length(barcodes.expBars)
%     curIdx = barcodes.idx(molIdx);
%     str = sprintf('Mol. %i',curIdx);
%     text(movies.pos{curIdx}(2),movies.pos{curIdx}(1)+voffList(molIdx),str,'Color','white','Clipping','on');
%  end

% figure
% for manual detection plot
[output,hPanelResult,images,movies,barcodes]  = sdd_script('sdd_fig_obed.txt');

molIdx = 19;

[ymin] = movies.pos{molIdx}(1);
[xmin] = movies.pos{molIdx}(2);
[ymax] = ymin+size(movies.molM{molIdx},1)-1;
[xmax] = xmin+size(movies.molM{molIdx},2)-1;

% nexttile
% imagesc(xmin:xmax,ymin:ymax,movies.molM{molIdx});colormap(gray)
% hold on
% plot(movies.trueedge{molIdx}(:, 2), movies.trueedge{molIdx}(:, 1),'red');
% 
% % plot(movies.trueedge{molIdx}(:, 2)-xmin+1, movies.trueedge{molIdx}(:, 1)-ymin+1,'red');
% title(['Mol.',num2str(molIdx)])

nexttile
imagesc(movies.dotM{molIdx});colormap(gray)
hold on
pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10)
axis off

% pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
% plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10)
for j = 1:length(barcodes.dots{molIdx}.val)
    str = sprintf('I = %.1f', barcodes.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
    text(barcodes.xy{molIdx}{2}(pos(j))-5,barcodes.xy{molIdx}{1}(pos(j))-5,str,'Color','white','Clipping','on');
end


hold on
ax2=plot(posx,posy,'green+','MarkerSize',10)

lgd = legend({'SDD detected dots','manual detection'})
lgd.Location ='southoutside';
% 
% f = figure
% imagesc(movies.dotM{molIdx});colormap(gray)
% [posx,posy] = getpts(f)

saveas(f,'Fig4.png');

