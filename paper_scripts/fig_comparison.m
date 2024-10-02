function [] = fig_comparison()

% data to analyze also with Tiff16

% first run SSD with settings from sdd_fig4.txt on one of the datasets
ix = 1;

datafolds = {'/export/scratch/albertas/data_temp/DOTS/data_for_figures/031523 h202 data/1/h202/','/export/scratch/albertas/data_temp/DOTS/data_for_figures/031523 h202 data/1/ut/','/export/scratch/albertas/data_temp/DOTS/data_for_figures/031523 h202 data/2/h202/','/export/scratch/albertas/data_temp/DOTS/data_for_figures/031523 h202 data/2/ut/' };
[output,hPanelResult,images,movies,barcodes]  = sdd_script('sdd_fig4.txt',[],datafolds{ix});

fold = output{1}.molRunFold;
files = dir(fullfile(fold,'*.tif')); % make sure the file ordering is the same!!
% fn = arrayfun(@(x) files(x).name,1:length(files),'un',false)';

fn = arrayfun(@(x) extract(files(x).name,digitsPattern),1:length(files),'un',false)';
mv=cellfun(@(x) str2num(x{1}),fn);
fn=cellfun(@(x) str2num(x{2}),fn);
% [a,b] = sortrows([mv fn]',2)';
[a,b] = sortrows([mv fn],'ascend');
% [a,b] = sort(fn);
% files = files(b);

% extract parameters from the results file
xlsFiles = dir([strrep(fold,'molecules_','results_'),'*.xlsx']);
calculatedLengths = zeros(1,length(files));
calculatedNumDots = zeros(1,length(files));
calculatedEdgeDots = zeros(1,length(files));
calculatedConvexToHull = zeros(1,length(files));

dataSSD = cell(1,length(xlsFiles));

ix = 1;
for i=1:length(xlsFiles)
    dataSSD{i} = readcell(fullfile(xlsFiles(i).folder,xlsFiles(i).name));
    calculatedLengths(ix:ix+size(dataSSD{i},1)-3) = cell2mat(dataSSD{i}(3:end,2));
        calculatedNumDots(ix:ix+size(dataSSD{i},1)-3) = cell2mat(dataSSD{i}(3:end,3));
    calculatedEdgeDots(ix:ix+size(dataSSD{i},1)-3) = cell2mat(dataSSD{i}(3:end,4));
    calculatedConvexToHull(ix:ix+size(dataSSD{i},1)-3) = cell2mat(dataSSD{i}(3:end,6));
    ix = ix+size(dataSSD{i},1)-2;
 end

%  calculatedLengths = calculatedLengths(b);
outname = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

foldRun1 = ['tifftest',outname];
mkdir(foldRun1);

cvThresh = 0.7; % save mols with specified parameters
molsKeep = calculatedConvexToHull>cvThresh;
ii=1;
for k=1:length(files)
    if molsKeep(k)==1

        A1 = imread(fullfile(files(k).folder,files(k).name),1);
        A2 = imread(fullfile(files(k).folder,files(k).name),2);
        name = strcat(['Green (', num2str(ii), ').tif']);
        imwrite(A1,fullfile(foldRun1,name));
        name = strcat(['Red (', num2str(ii), ').tif']);
        imwrite(A2,fullfile(foldRun1,name));
        ii=ii+1;
    end
end


foldRun = ['sddsmall',outname];
mkdir(foldRun);
ii=1;
for k=1:length(files)
    if molsKeep(k)==1
        A1 = imread(fullfile(files(k).folder,files(k).name),1);
        A2 = imread(fullfile(files(k).folder,files(k).name),2);
        name = strcat(['Green (', num2str(ii), ').tiff']);
        imwrite(A1,fullfile(foldRun,name));
        imwrite(A2,fullfile(foldRun,name),'WriteMode', "append");
        ii=ii+1;
    end
end

%% Now run the Tiff16 analizer
foldRun1 = '/export/scratch/albertas/data_temp/DOTS/oUTPUT/ix1/tifftest2024-05-15_15_42_43_test1_output_ix1/tifftest2024-05-15_15_42_43/';
foldRun= '/export/scratch/albertas/data_temp/DOTS/code_for_figures/sddsmall2024-05-15_15_42_43/';

% import results from Tiff16 analizer
A = importdata(fullfile(foldRun1,'result.txt'));

infoLine = find(cellfun(@(x) contains(x,'ImageNum'),A));
infoText = A(infoLine);

data = A(infoLine+2:end);
% data(output{1}.delid) = [];

mols.imageNum =1:length(data);
mols.greenLength = zeros(1,length(data));
for i=1:length(data)
    splitData = strsplit(data{i}, '|');
    mols.greenLength(i) = str2num(splitData{2});
    mols.redNum(i) = str2num(splitData{4});
end

%%
mols.greenLength = mols.greenLength(b);
mols.redNum = mols.redNum(b);
keepRows =  find(mols.greenLength>=40);
lengthsM =mols.greenLength*110/1000;

tiff16_lengths = lengthsM(keepRows);
sdd_lengths = calculatedLengths(keepRows);

lenDigLessThanMicron = find(abs(tiff16_lengths-sdd_lengths)<1);
% calculatedLengths(keepRows),lengthsM(keepRows)


sddDots = calculatedNumDots(keepRows);
tiff16Dots =  mols.redNum(keepRows);

sddDots = sddDots(lenDigLessThanMicron);
tiff16Dots = tiff16Dots(lenDigLessThanMicron);

keptMols = keepRows(lenDigLessThanMicron);
[sddDots;tiff16Dots]
% keptMols(139)
%Engine
[uxy, jnk, idx] = unique([sddDots.',tiff16Dots.'],'rows');
szscale = histcounts(idx,min(unique(idx))-0.5:max(unique(idx))+0.5);
%Plot Scale of 25 and stars
f=figure
hold on
ax = scatter(uxy(:,1),uxy(:,2),'c','blue','filled','sizedata',szscale*3)
ax.MarkerEdgeColor = [0.85  0.32 0.098];
ax.MarkerFaceColor = [0.85  0.32 0.098];
xlabel('SDD Number of Dots')
ylabel('Tiff16 number of dots')

mdl = fitlm(sddDots,tiff16Dots,'intercept',false)
plot(min(sddDots):max(sddDots),mdl.Coefficients.Estimate(1)*(min(sddDots):max(sddDots)))
lgd=legend({['$R^2$ =',num2str(mdl.Rsquared.Ordinary)],['f(x) = ',num2str(mdl.Coefficients.Estimate(1))]},'Interpreter','latex')
print(f,['FigTiff16dots_' num2str(ix) '.png'], '-dpng', '-r300', '-painters');

% plot(min(x):max(x),mdl.Coefficients.Estimate(2)*(min(x):max(x))+mdl.Coefficients.Estimate(1))
% lgd=legend({['$R^2$ =',num2str(mdl.Rsquared.Ordinary)],['f(x) = ',num2str(mdl.Coefficients.Estimate(2)),'x','+',num2str(mdl.Coefficients.Estimate(1))]},'Interpreter','latex')
% lgd.Location ='eastoutside';

%%
% results from both. Tiff16 - import, sdd - detect (for this molecule
% specifically)
% idx = b(end); % good
% goodMols = b(keepRows);
idx=143
[B, out1,hPanelResult1,images1,movies1,barcodes1] = tiff16_vs_sdd_plot(foldRun1,foldRun,b,keepRows,idx);

idx = 7; % good
[B, out1,hPanelResult1,images1,movies1,barcodes1] = tiff16_vs_sdd_plot(foldRun1,foldRun,b,keepRows,idx);

% 
% % figure;imshow(B)
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


%%
f=figure
tiledlayout(3, 2,'TileSpacing','tight','Padding','none')
ax1 = nexttile
hold on
title('A) Tiff16analizer output','Interpreter','latex','FontName','Times')
imagesc(B)
axis off

axis equal;            % <---- move to after-plot
daspect(ax1,[1 1 1]);  % <---- move to after-plot
pbaspect(ax1,[1 0.32 1]); % <---- move to after-plot
% 
xbound = size(B,2);
ybound = size(B,1);

nPixels = 1e3/sets.pxnm; % 1 micron
x = [10 nPixels+10];
y = [3 3 ];
plot(x,y,'Linewidth',2,'Color','white')

ax2=nexttile

% show all
hold on
title('B) SDD output','Interpreter','latex','FontName','Times')
imagesc(images1{1}.dotIm);colormap(gray)
hold on

xbound = size(images1{1}.dotIm,2);
ybound = size(images1{1}.dotIm,1);

nPixels = 1e3/sets.pxnm;
x = [10 nPixels+10];
y = [3 3 ];
plot(x,y,'Linewidth',2,'Color','white')


axis equal;            % <---- move to after-plot
daspect(ax2,[1 1 1]);  % <---- move to after-plot
pbaspect(ax2,[1 0.32 1]); % <---- move to after-plot
axis off

for molIdx=1:length(movies1.dotM)
    plot(movies1.trueedge{molIdx}(:, 2), movies1.trueedge{molIdx}(:, 1),'green');
    pos = barcodes1.dots{molIdx}.locations+barcodes1.dots{molIdx}.leftOffset;

    plot(barcodes1.xy{molIdx}{2}(pos),barcodes1.xy{molIdx}{1}(pos),'redx','MarkerSize',10)
    for j = 1:length(barcodes1.dots{molIdx}.val)
        str = sprintf('I = %.1f', barcodes1.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
        text(barcodes1.xy{molIdx}{2}(pos(j))-5,barcodes1.xy{molIdx}{1}(pos(j))-5,str,'Color','white','Clipping','on','FontSize',8);
    end
end

% axis off
% imagesc(movies1.dotM{molIdx});colormap(gray)
% hold on
% pos = barcodes1.dots{molIdx}.locations+barcodes1.dots{molIdx}.leftOffset;
% %     pos
% plot(barcodes1.xy{molIdx}{2}(pos),barcodes1.xy{molIdx}{1}(pos),'redx','MarkerSize',10)
% for j = 1:length(barcodes1.dots{molIdx}.val)
%     str = sprintf('I = %.1f', barcodes1.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
%     text(barcodes1.xy{molIdx}{2}(pos(j))-5,barcodes1.xy{molIdx}{1}(pos(j))-5,str,'Color','white','Clipping','on');
% end
ax3 = nexttile([2 1])

scatter(calculatedLengths(keepRows),lengthsM(keepRows))
xlabel('Length SDD ($\mu m$)','Interpreter','latex','FontName','Times')
ylabel('Length Tiff16 ($\mu m$)','Interpreter','latex','FontName','Times')
title('C) Lengths SDD vs Tiff16','Interpreter','latex','FontName','Times')
hold on
axis equal
% daspect(ax3,[1 1 1]);  % <---- move to after-plot
% pbaspect(ax3,[1 1 1]); % <---- move to after-plot
% daspect(ax1,[1 1 1]);  % <---- move to after-plot
% pbaspect([2.4 0.5 1])

y = lengthsM(keepRows);
x = calculatedLengths(keepRows);

mdl = fitlm(x,y)
% text(1,max(y)-1,['f(x) = ',num2str(mdl.Coefficients.Estimate(2)),'x','+',num2str(mdl.Coefficients.Estimate(1))])
% text(1,max(y)-2,['R^2 =',num2str(mdl.Rsquared.Ordinary)])

plot(min(x):max(x),mdl.Coefficients.Estimate(2)*(min(x):max(x))+mdl.Coefficients.Estimate(1))
if mdl.Coefficients.Estimate(1) < 1
    lgd=legend({['$R^2$ =',num2str(mdl.Rsquared.Ordinary)],['f(x) = ',num2str(mdl.Coefficients.Estimate(2)),'x',num2str(mdl.Coefficients.Estimate(1))]},'Interpreter','latex')
else
    lgd=legend({['$R^2$ =',num2str(mdl.Rsquared.Ordinary)],['f(x) = ',num2str(mdl.Coefficients.Estimate(2)),'x','+',num2str(mdl.Coefficients.Estimate(1))]},'Interpreter','latex')
end
lgd.Location ='southoutside';

[y(idx) x(idx)]
% nexttile
% molIdx = 1;
% imagesc(movies1.molM{molIdx});colormap(gray)
% hold on
% plot(movies1.trueedge{molIdx}(:, 2)-xmin+1, movies1.trueedge{molIdx}(:, 1)-ymin+1,'red');
ax1=nexttile([2 1])
title('D) Dots SDD vs Tiff16','Interpreter','latex','FontName','Times')
%Engine
[uxy, jnk, idx] = unique([sddDots.',tiff16Dots.'],'rows');
szscale = histcounts(idx,min(unique(idx))-0.5:max(unique(idx))+0.5);
%Plot Scale of 25 and stars
% f=figure
hold on
ax = scatter(uxy(:,1),uxy(:,2),'c','blue','filled','sizedata',szscale*3)
ax.MarkerEdgeColor = [0.85  0.32 0.098];
ax.MarkerFaceColor = [0.85  0.32 0.098];
xlabel('SDD Number of Dots','Interpreter','latex','FontName','Times')
ylabel('Tiff16 number of dots','Interpreter','latex','FontName','Times')
axis equal;            % <---- move to after-plot

xlim([0 max(uxy(:,1)) ])
ylim([0 max(uxy(:,2)) ])

mdl = fitlm(tiff16Dots,sddDots,'intercept',false)
plot(min(sddDots):max(sddDots),mdl.Coefficients.Estimate(1)*(min(sddDots):max(sddDots)))
lgd=legend({['$R^2$ =',num2str(mdl.Rsquared.Ordinary)],['f(x) = ',num2str(mdl.Coefficients.Estimate(1)),'x']},'Interpreter','latex')
% print(f,['FigTiff16dots_' num2str(ix) '.png'], '-dpng', '-r300', '-painters');
lgd.Location ='southoutside';

% daspect(ax1,[1 1 1]);  % <---- move to after-plot
% pbaspect(ax1,[1 1 1]); % <---- move to after-plot
set(gcf, 'Color', 'w')

print(f, 'FigTiff16.png', '-dpng', '-r300', '-painters');
%%
figure,
hold on
bar(1,mean(calculatedLengths(keepRows)))                
errorbar(mean(calculatedLengths(keepRows)), std(calculatedLengths(keepRows)))


%%
foldRun = ['tiffDataTest',outname];
mkdir(foldRun);
ii=1;
for k=1:length(keepRows)
        A1 = imread(fullfile(files(keepRows(k)).folder,files(keepRows(k)).name),1);
        A2 = imread(fullfile(files(keepRows(k)).folder,files(keepRows(k)).name),2);
        name = strcat(['Mol (', num2str(ii), ').tiff']);
        imwrite(A1,fullfile(foldRun,name));
        imwrite(A2,fullfile(foldRun,name),'WriteMode', "append");
        ii=ii+1;
   
end



end

