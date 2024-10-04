function [] = fig_dots()

[output,hPanelResult,images,movies,barcodes]  = sdd_script('sdd_fig_obed.txt');

%%

f=figure
tiledlayout(3,3,'TileSpacing','compact');
nexttile([3 2])
imagesc(images{1}.registeredIm{1});colormap(gray)

molIdx = 18;

[ymin] = movies.pos{molIdx}(1);
[xmin] = movies.pos{molIdx}(2);
[ymax] = ymin+size(movies.molM{molIdx},1)-1;
[xmax] = xmin+size(movies.molM{molIdx},2)-1;


rectangle('Position', [xmin, ymin, size(movies.molM{molIdx},2), size(movies.molM{molIdx},1)], 'EdgeColor', 'r', 'LineWidth', 1);

% imshowpair(images{1}.registeredIm{1},images{1}.dotIm);colormap(gray)
% axis off
hold on
% line([xmin, size(images{1}.registeredIm{1},2)], [ymin+size(movies.molM{molIdx},1), size(images{1}.registeredIm{1},1)/3], 'Color', 'r', 'LineWidth', 1);
% line([xmin, size(images{1}.registeredIm{1},2)], [ymin+size(movies.molM{molIdx},1), 1], 'Color', 'r', 'LineWidth', 1);

nexttile
imagesc(movies.molM{molIdx});colormap(gray)
hold on
plot(movies.trueedge{molIdx}(:, 2)-xmin+1, movies.trueedge{molIdx}(:, 1)-ymin+1,'red');

nexttile
imagesc(movies.dotM{molIdx});colormap(gray)
hold on
pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10)
pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10)
for j = 1:length(barcodes.dots{molIdx}.val)
    str = sprintf('I = %.1f', barcodes.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
    text(barcodes.xy{molIdx}{2}(pos(j))-5,barcodes.xy{molIdx}{1}(pos(j))-5,str,'Color','white','Clipping','on');
end
print('FigDots.eps','-depsc','-r300');
saveas(f,'FigDots.png');

 %%

hold on
trueedge = movies.trueedge;
%     axes(tiles.dotDet);
hold on

for k = 1:length(trueedge)
%       figure(dotFigNum)
  plot(trueedge{k}(:, 2), trueedge{k}(:, 1));
end

for molIdx = 1:length(barcodes.dots)
    curIdx = barcodes.idx(molIdx);
    str = sprintf('Mol. %i',curIdx);
    angle = atan(barcodes.lineParams{curIdx}(1));
    % 
    vOff = barcodes.boundaries{curIdx}(1);
    hOff = barcodes.boundaries{curIdx}(3);
    voffList(molIdx) = vOff;
end

hold on
 for molIdx = 1:length(barcodes.expBars)
    curIdx = barcodes.idx(molIdx);
    str = sprintf('Mol. %i',curIdx);
    text(movies.pos{curIdx}(2),movies.pos{curIdx}(1)+voffList(molIdx),str,'Color','white','Clipping','on');
 end

 mIdx = [14 19 11]

 cnvToHull = cellfun(@(x) x.FilledArea/x.ConvexArea,movies.stats);

 ecc = cellfun(@(x) x.Eccentricity,movies.stats);


lengths = cellfun(@(x) length(x.rawBarcode),barcodes.expBars);



rat = 0:0.01:1;
 for i= 1:length(rat)
     numRem(i) = sum(cnvToHull>rat(i));
    eccrem(i) = sum(ecc>rat(i));
    lengthsrem(i) = sum(lengths > rat(i)*max(lengths));
 end

xlim([600 1100])
ylim([100 800])

nPixels = 1e4/sets.pxnm;
x = 600+[5, 5 + nPixels ];
y = [0.9*800 ,0.9*800 ];
plot(x,y,'Linewidth',8,'Color','white')
text(0,0.05,'10 microns','Fontsize',10,'Color','white','Units','normalized')

axis off

nexttile
hold on
plot(rat,numRem/length(cnvToHull))
title('Convex to hull ratio')
xlabel('ratlim')
ix = 26;
plot(cnvToHull(ix),numRem(find(cnvToHull(ix)<rat,1,'first'))/length(cnvToHull),'redx')

% ylabel('Remaining barcodes')



nexttile
hold on
plot(rat,eccrem/length(cnvToHull))
title('Eccentricity')
xlabel('elim')
% ylabel('Remaining barcodes')
plot(ecc(ix),eccrem(find(ecc(ix)<rat,1,'first'))/length(cnvToHull),'redx')




nexttile
hold on
plot( rat*max(lengths),lengthsrem/length(cnvToHull))
title('Length (px) ')
xlabel('Length (px)')
ylabel('Remaining barcodes')

plot(lengths(ix),lengthsrem(find(lengths(ix)<rat*max(lengths),1,'first'))/length(cnvToHull),'redx')

%%
sets.pxnm = output{1}.settings.pixelSize;

f=figure
% tiledlayout(4,3,'TileSpacing','compact');
% nexttile([4 2])
imagesc(images{1}.dotIm);colormap(gray)
hold on
    axis equal


trueedge = movies.trueedge;
%     axes(tiles.dotDet);
hold on



for k = 1:length(trueedge)
%       figure(dotFigNum)
  plot(trueedge{k}(:, 2), trueedge{k}(:, 1));
end


molIdx = 13;

[ymin] = movies.pos{molIdx}(1);
[xmin] = movies.pos{molIdx}(2);
[ymax] = ymin+size(movies.molM{molIdx},1)-1;
[xmax] = xmin+size(movies.molM{molIdx},2)-1;
rectangle('Position', [xmin, ymin, size(movies.molM{molIdx},2), size(movies.molM{molIdx},1)], 'EdgeColor', 'cyan', 'LineWidth', 1);



molIdx = 11;
[ymin] = movies.pos{molIdx}(1);
[xmin] = movies.pos{molIdx}(2);
[ymax] = ymin+size(movies.molM{molIdx},1)-1;
[xmax] = xmin+size(movies.molM{molIdx},2)-1;
% rectangle('Position', [xmin, ymin, size(movies.molM{molIdx},2), size(movies.molM{molIdx},1)], 'EdgeColor', 'r', 'LineWidth', 1);

% 
sizeY1 = max(size(movies.dotM{17},1),size(movies.dotM{19},1));
sizeX1 = max(size(movies.dotM{19},2),size(movies.dotM{17},2));

for molIdxT = 1:length(barcodes.dots)
    curIdx = barcodes.idx(molIdxT);
    str = sprintf('Mol. %i',curIdx);
    angle = atan(barcodes.lineParams{curIdx}(1));
    % 
    vOff = barcodes.boundaries{curIdx}(1);
    hOff = barcodes.boundaries{curIdx}(3);
    voffList(molIdxT) = vOff;
end

hold on

 for molIdxT = 1:length(barcodes.expBars)
    curIdx = barcodes.idx(molIdxT);
    str = sprintf('Mol. %i',curIdx);
    text(movies.pos{curIdx}(2),movies.pos{curIdx}(1)+voffList(molIdxT),str,'Color','white','Clipping','on');
 end

axis off
% xbound = size(images{1}.dotIm,2);
% ybound =  size(images{1}.dotIm,1);
% nPixels = 1e3/sets.pxnm;
% x = [xbound-10*nPixels-10 xbound-10];
% y = [ybound-20 ybound-20 ];
% plot(x,y,'Linewidth',2,'Color','white')
xlim([430 966])
ylim([394 796])

xbound = 966;
ybound =  796;
nPixels = 1e3/sets.pxnm;
x = [xbound-10*nPixels-10 xbound-10];
y = [ybound-20 ybound-20 ];
plot(x,y,'Linewidth',2,'Color','white')


set(gcf, 'Color', 'w')

print('FigDotsD1.png','-dpng','-r400');
%
f=figure
tiledlayout(2,2,'TileSpacing','compact');

nexttile
imagesc(xmin:xmax,ymin:ymax,movies.molM{molIdx});colormap(gray)
hold on
plot(movies.trueedge{molIdx}(:, 2), movies.trueedge{molIdx}(:, 1),'red');
    axis equal

% plot(movies.trueedge{molIdx}(:, 2)-xmin+1, movies.trueedge{molIdx}(:, 1)-ymin+1,'red');
title(['Mol.',num2str(molIdx)])
axis off
x = [xmin+10 xmin+nPixels+10];
y = [ymax-5 ymax-5 ];
plot(x,y,'Linewidth',2,'Color','white')


nexttile
imagesc(movies.dotM{molIdx});colormap(gray)
hold on
    axis equal

pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10,'LineWidth',1.5)
axis off

pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10,'LineWidth',1.5)
% for j = 1:length(barcodes.dots{molIdx}.val)
%     str = sprintf('I = %.1f', barcodes.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
%     text(barcodes.xy{molIdx}{2}(pos(j))-5,barcodes.xy{molIdx}{1}(pos(j))-5,str,'Color','white','Clipping','on');
% end

% Another molecule
molIdx = 13;
[ymin] = movies.pos{molIdx}(1);
[xmin] = movies.pos{molIdx}(2);
[ymax] = ymin+size(movies.molM{molIdx},1)-1;
[xmax] = xmin+size(movies.molM{molIdx},2)-1;
nexttile
imagesc(xmin:xmax,ymin:ymax,movies.molM{molIdx});colormap(gray)
hold on
plot(movies.trueedge{molIdx}(:, 2), movies.trueedge{molIdx}(:, 1),'cyan');
axis equal

% plot(movies.trueedge{molIdx}(:, 2)-xmin+1, movies.trueedge{molIdx}(:, 1)-ymin+1,'red');
title(['Mol.',num2str(molIdx)])
axis off
x = [xmin+10 xmin+nPixels+10];
y = [ymax-2 ymax-2 ];
plot(x,y,'Linewidth',2,'Color','white')


nexttile
imagesc(movies.dotM{molIdx});colormap(gray)
hold on
    axis equal

pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'cyanx','MarkerSize',10,'LineWidth',1.5)
axis off

pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'cyanx','MarkerSize',10,'LineWidth',1.5)
% for j = 1:length(barcodes.dots{molIdx}.val)
%     str = sprintf('I = %.1f', barcodes.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
%     text(barcodes.xy{molIdx}{2}(pos(j))-5,barcodes.xy{molIdx}{1}(pos(j))-5,str,'Color','white','Clipping','on');
% end

set(gcf, 'Color', 'w')

% ylim([1 sizeY1])
%



print('FigDotsD2.png','-dpng','-r400');

%% SEP


%%
sets.pxnm = output{1}.settings.pixelSize;


molIdx = 17;
[ymin] = movies.pos{molIdx}(1);
[xmin] = movies.pos{molIdx}(2);
[ymax] = ymin+size(movies.molM{molIdx},1)-1;
[xmax] = xmin+size(movies.molM{molIdx},2)-1;
% rectangle('Position', [xmin, ymin, size(movies.molM{molIdx},2), size(movies.molM{molIdx},1)], 'EdgeColor', 'r', 'LineWidth', 1);

% 
sizeY1 = max(size(movies.dotM{17},1),size(movies.dotM{19},1));
sizeX1 = max(size(movies.dotM{19},2),size(movies.dotM{17},2));

for molIdxT = 1:length(barcodes.dots)
    curIdx = barcodes.idx(molIdxT);
    str = sprintf('Mol. %i',curIdx);
    angle = atan(barcodes.lineParams{curIdx}(1));
    % 
    vOff = barcodes.boundaries{curIdx}(1);
    hOff = barcodes.boundaries{curIdx}(3);
    voffList(molIdxT) = vOff;
end

% hold on

%  for molIdxT = 1:length(barcodes.expBars)
%     curIdx = barcodes.idx(molIdxT);
%     str = sprintf('Mol. %i',curIdx);
%     text(movies.pos{curIdx}(2),movies.pos{curIdx}(1)+voffList(molIdxT),str,'Color','white','Clipping','on');
%  end
% 
% axis off
xbound = size(images{1}.dotIm,2);
ybound =  size(images{1}.dotIm,1);
nPixels = 1e3/sets.pxnm;
x = [xbound-10*nPixels-10 xbound-10];
y = [ybound-20 ybound-20 ];
% plot(x,y,'Linewidth',2,'Color','white')

figure
imagesc(xmin:xmax,ymin:ymax,movies.molM{molIdx});colormap(gray)
hold on
plot(movies.trueedge{molIdx}(:, 2), movies.trueedge{molIdx}(:, 1),'red');
    axis equal

% plot(movies.trueedge{molIdx}(:, 2)-xmin+1, movies.trueedge{molIdx}(:, 1)-ymin+1,'red');
title(['Mol.',num2str(molIdx)])
axis off
x = [xmin+10 xmin+nPixels+10];
y = [ymax-5 ymax-5 ];
plot(x,y,'Linewidth',2,'Color','white')

print('FigDots1D1.png','-dpng','-r400');

figure
imagesc(movies.dotM{molIdx});colormap(gray)
hold on
    axis equal

pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10)
axis off

pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10)
% for j = 1:length(barcodes.dots{molIdx}.val)
%     str = sprintf('I = %.1f', barcodes.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
%     text(barcodes.xy{molIdx}{2}(pos(j))-5,barcodes.xy{molIdx}{1}(pos(j))-5,str,'Color','white','Clipping','on');
% end
print('FigDots1D2.png','-dpng','-r400');

% Another molecule
molIdx = 19;
[ymin] = movies.pos{molIdx}(1);
[xmin] = movies.pos{molIdx}(2);
[ymax] = ymin+size(movies.molM{molIdx},1)-1;
[xmax] = xmin+size(movies.molM{molIdx},2)-1;
figure
imagesc(xmin:xmax,ymin:ymax,movies.molM{molIdx});colormap(gray)
hold on
plot(movies.trueedge{molIdx}(:, 2), movies.trueedge{molIdx}(:, 1),'cyan');
axis equal

% plot(movies.trueedge{molIdx}(:, 2)-xmin+1, movies.trueedge{molIdx}(:, 1)-ymin+1,'red');
title(['Mol.',num2str(molIdx)])
axis off
x = [xmin+10 xmin+nPixels+10];
y = [ymax-2 ymax-2 ];
plot(x,y,'Linewidth',2,'Color','white')

print('FigDots1D3.png','-dpng','-r400');

figure
imagesc(movies.dotM{molIdx});colormap(gray)
hold on
    axis equal

pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'cyanx','MarkerSize',10)
axis off

pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'cyanx','MarkerSize',10)
% for j = 1:length(barcodes.dots{molIdx}.val)
%     str = sprintf('I = %.1f', barcodes.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
%     text(barcodes.xy{molIdx}{2}(pos(j))-5,barcodes.xy{molIdx}{1}(pos(j))-5,str,'Color','white','Clipping','on');
% end

set(gcf, 'Color', 'w')
print('FigDots1D4.png','-dpng','-r400');

% ylim([1 sizeY1])
%





% ylim([1 sizeY1])

% print('FigDots.eps','-depsc','-r300');
% saveas(f,'panelC.png');
end

