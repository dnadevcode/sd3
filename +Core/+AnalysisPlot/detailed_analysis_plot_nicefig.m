function [] = detailed_analysis_plot_nicefig(images, movies,barcodes,idx)



f=figure;
tiledlayout(2,3,'TileSpacing','compact');
nexttile([2 2])
imagesc(images{1}.dotIm);colormap(gray)
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


% molIdx = 19;
molIdx = idx;
[ymin] = movies.pos{molIdx}(1);
[xmin] = movies.pos{molIdx}(2);
[ymax] = ymin+size(movies.molM{molIdx},1)-1;
[xmax] = xmin+size(movies.molM{molIdx},2)-1;

nexttile
imagesc(xmin:xmax,ymin:ymax,movies.molM{molIdx});colormap(gray)
hold on
plot(movies.trueedge{molIdx}(:, 2), movies.trueedge{molIdx}(:, 1),'red');

% plot(movies.trueedge{molIdx}(:, 2)-xmin+1, movies.trueedge{molIdx}(:, 1)-ymin+1,'red');
title(['Mol.',num2str(molIdx)])

nexttile
imagesc(movies.dotM{molIdx});colormap(gray)
hold on
pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10);
axis off

pos = barcodes.dots{molIdx}.locations+barcodes.dots{molIdx}.leftOffset;
%     pos
plot(barcodes.xy{molIdx}{2}(pos),barcodes.xy{molIdx}{1}(pos),'redx','MarkerSize',10);
for j = 1:length(barcodes.dots{molIdx}.val)
    str = sprintf('I = %.1f', barcodes.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
    text(barcodes.xy{molIdx}{2}(pos(j))-5,barcodes.xy{molIdx}{1}(pos(j))-5,str,'Color','white','Clipping','on');
end
