function [f] = fig3(images,movies,barcodes,sets)

lw = 2; % linewidth
f=figure
tiledlayout(3,3,'TileSpacing','compact');
nexttile([3 2])

imagesc(images{1}.registeredIm{1});colormap(gray)
% axis off
hold on


hold on
trueedge = movies.trueedge;
%     axes(tiles.dotDet);
hold on

ix1 = 26;
ix2 = 28;

for k = 1:length(trueedge)
%       figure(dotFigNum)
 if k==ix1
    plot(trueedge{k}(:, 2), trueedge{k}(:, 1),'red');
 else 
     if k==ix2
    plot(trueedge{k}(:, 2), trueedge{k}(:, 1),'blue');
 else
    plot(trueedge{k}(:, 2), trueedge{k}(:, 1),'green');
     end
 end
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
    text(movies.pos{curIdx}(2),movies.pos{curIdx}(1)+voffList(molIdx),str,'Fontsize',8,'Color','white','Clipping','on');
 end

%  mIdx = [14 19 11]
 cnvToHull = cellfun(@(x) x.FilledArea/x.ConvexArea,movies.stats);

 ecc = cellfun(@(x) x.Eccentricity,movies.stats);

 lengths = nan(1,length(movies.pos));
allBars = 1:length(movies.pos);
allBars(barcodes.delid) = 0;
lengths(logical(allBars)) = cellfun(@(x) length(x.rawBarcode),barcodes.expBars);



rat = 0:0.01:1;
 for i= 1:length(rat)
     numRem(i) = sum(cnvToHull>rat(i));
    eccrem(i) = sum(ecc>rat(i));
    lengthsrem(i) = sum(lengths > rat(i)*max(lengths));
 end

 xbound = [600 1100];
xlim(xbound);
ybound = [100 800];
ylim(ybound);

nPixels = 1e4/sets.pxnm;
x = [xbound(2)-nPixels-10 xbound(2)-10];
y = [ybound(2)-20 ybound(2)-20 ];
plot(x,y,'Linewidth',2,'Color','white')
% text(0,0.05,'10 microns','Fontsize',10,'Color','white','Units','normalized')
% title('A) SDD output')
axis off

nexttile
hold on
plot(rat,numRem/length(cnvToHull),'black')
% title('B) Convex to hull ratio')
xlabel('$ratlim_{thresh}$','Interpreter','latex')


plot(cnvToHull(ix1),numRem(find(cnvToHull(ix1)<rat,1,'first'))/length(cnvToHull),'redx','LineWidth',lw)
plot(cnvToHull(ix2),numRem(find(cnvToHull(ix2)<rat,1,'first'))/length(cnvToHull),'bluex','LineWidth',lw)

plot([0.4 0.4],[0 1], 'red-','LineWidth',lw)
% ylabel('Remaining barcodes')



nexttile
hold on
plot(rat,eccrem/length(cnvToHull),'black')
% title('C) Eccentricity')
xlabel('$elim_{thresh}$','Interpreter','latex')
% ylabel('Remaining barcodes')
plot(ecc(ix1),eccrem(find(ecc(ix1)<rat,1,'first'))/length(cnvToHull),'redx','LineWidth',lw)
plot(ecc(ix2),eccrem(find(ecc(ix2)<rat,1,'first'))/length(cnvToHull),'bluex','LineWidth',lw)
plot([0.8 0.8],[0 1], 'red-','LineWidth',lw)




nexttile
hold on
plot(lengths(ix1),lengthsrem(find(lengths(ix1)<rat*max(lengths),1,'first'))/length(cnvToHull),'redx','LineWidth',lw)
plot(lengths(ix2),lengthsrem(find(lengths(ix2)<rat*max(lengths),1,'first'))/length(cnvToHull),'bluex','LineWidth',lw)
plot([50 50],[0 1], 'red-','LineWidth',lw)

plot( rat*max(lengths),lengthsrem/length(cnvToHull),'black')

% title('D) Length (px) ')
xlabel('$Length_{thresh}$','Interpreter','latex')
ylabel('Fraction of remaining molecules','FontSize',10,'FontName','Times')


lgd=legend({['Mol.',num2str(ix1)],['Mol.',num2str(ix2)],'Default threshold'})
lgd.Location='southoutside';