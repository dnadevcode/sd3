% mols prep for tiffanalizer

%{
 cd "/export/scratch/albertas/data_temp/DOTS/031523 h202 data/1/h202/molecules_run3"
lcd C:\Users\Lenovo\postdoc\DATA\DOT\OneDrive_2023-10-19\out
%}

fold = 'C:\Users\Lenovo\postdoc\DATA\DOT\OneDrive_2023-10-19\out\';


files = dir(fullfile(fold,'*.tif'));

foldRun = 'tifftest';
mkdir(foldRun);


for i=1:length(files)

    A1 = imread(fullfile(files(i).folder,files(i).name),1);
    A2 = imread(fullfile(files(i).folder,files(i).name),2);

    name = strcat(['Green (', num2str(i), ').tif']);
    imwrite(A1,fullfile(foldRun,name));
    name = strcat(['Red (', num2str(i), ').tif']);
    imwrite(A2,fullfile(foldRun,name));
end




xlsFiles = dir(fullfile(fileparts(fold),'\*.xlsx'));
calculatedLengths = [];
calculatedNumDots = [];
calculatedEdgeDots = [];

for i=1:length(xlsFiles)
    dataSSD{i} = readcell(fullfile(xlsFiles{i}.folder,xlsFiles{i}.name));
    calculatedLengths = [calculatedLengths cell2mat(dataSSD{i}(3:end,2))];
    calculatedNumDots = [calculatedNumDots cell2mat(dataSSD{i}(3:end,3))];
    calculatedEdgeDots= [calculatedEdgeDots cell2mat(dataSSD{i}(3:end,4))];
end





load('C:\Users\Lenovo\postdoc\DATA\DOT\OneDrive_2023-10-19\out\dnarecoutput_1.mat')
A = importdata(fullfile(foldRun,'result.txt'));

infoLine = find(cellfun(@(x) contains(x,'ImageNum'),A));
infoText = A(infoLine);

data = A(infoLine+2:end);
data(output{1}.delid) = [];

mols.imageNum =1:length(data);
mols.greenLength = zeros(1,length(data));
for i=1:length(data)
    splitData = strsplit(data{i}, '|');
    mols.greenLength(i) = str2num(splitData{2}) 
    mols.redNum(i) = str2num(splitData{4}) 
end


xlsFiles = dir('C:\Users\Lenovo\postdoc\DATA\DOT\OneDrive_2023-10-19\out\*.xlsx');
calculatedLengths = [];
calculatedNumDots = [];
calculatedEdgeDots = [];

for i=1:length(xlsFiles)
    dataSSD{i} = readcell(fullfile(xlsFiles{i}.folder,xlsFiles{i}.name));
    calculatedLengths = [calculatedLengths cell2mat(dataSSD{i}(3:end,2))];
    calculatedNumDots = [calculatedNumDots cell2mat(dataSSD{i}(3:end,3))];
    calculatedEdgeDots= [calculatedEdgeDots cell2mat(dataSSD{i}(3:end,4))];
end




keepRows =  find(mols.greenLength~=0);
figure,plot( mols.redNum(keepRows),calculatedNumDots(keepRows),'x')

x = calculatedNumDots(keepRows)';
y =  mols.redNum(keepRows);
%Engine
[uxy, jnk, idx] = unique([x.',y.'],'rows');
szscale = histcounts(idx,min(unique(idx))-0.5:max(unique(idx))+0.5);
%Plot Scale of 25 and stars
figure
scatter(uxy(:,1),uxy(:,2),'c','red','filled','sizedata',szscale*25)


B = imread('C:\Users\Lenovo\postdoc\DATA\DOT\OneDrive_2023-10-19\tifftest\RG (145).jpeg');
i=145;
 A1 = imread(fullfile(files(i).folder,files(i).name),1);


 lengthsM =mols.greenLength*110/1000;
% figure;imshow(B)
figure
tiledlayout(2,2)
nexttile
imagesc(A1)
nexttile
imagesc(B);colormap(gray)
title('Mol 145')
nexttile
scatter(lengthsM(keepRows),calculatedLengths(keepRows))
xlabel('Length SDD')
ylabel('Length tiff16')
nexttile
ax=scatter(uxy(:,1),uxy(:,2),'filled','c','sizedata',szscale*25)
ax.MarkerEdgeColor = [0.85  0.32 0.098];
ax.MarkerFaceColor = [0.85  0.32 0.098];
xlabel('Num dots SDD')
ylabel('Num dots tiff16')


    

