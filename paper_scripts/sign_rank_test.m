

xlsFiles = dir(['/export/scratch/albertas/data_temp/DOTS/data_for_figures/ManualLength/*.xlsx']);
excelFilePath = fullfile(xlsFiles(1).folder,xlsFiles(1).name);

sheetNamesXLS = sheetnames(excelFilePath);

% sheetNames = xlsfinfo(excelFilePath);
% dataSSDmanual = readtable(excelFilePath);

dataSSDmanual = cell(1,length(sheetNamesXLS))

for i=1:length(sheetNamesXLS)

    dataSSDmanual{i} = readcell(excelFilePath,'Sheet',sheetNamesXLS{i});

    x = cell2mat(dataSSDmanual{i}(2:end,1));
    y = cell2mat(dataSSDmanual{i}(2:end,2))
    mdl{i} = fitlm(x,y)
    
        p{i} = signrank(x,y)
end
%     dotsSDD = cell2mat(dataSSDmanual{i}(3:end,1));
%     dotsManual = cell2mat(dataSSDmanual{i}(3:end,2));
% end


