function [outputNew, allKymos] = goodbadtool(numImages,fold, data, foldOut)
    %   Args:
    %   numImages - array of number of images to display in x and y ,
    %   fold - folder with input images
    %   foldOut - output for folders good and bad molecules,foldOut
    %
    %   Saves all good files into good folder, bad into bad folder
    
    if nargin <2 
        fold = 'testfolder\molecules_run1';
        data = 'testfolder\dnarecoutput';
        foldOut = 'output';
    end
    
    load(data)
    import UI.select_image;
    
    if nargin< 1
        numImages = [5 1]; % grid for images
    end

    if nargin < 2
        % folder with images to classify
        fold = uigetdir(pwd, "Folder with tifs we want to classify");
    end


    if nargin < 3
        % folder with images to classify
        foldOut = uigetdir(pwd,"Select output folder");
    end

    numImagesToShow = numImages(1)*numImages(2);

    listing = [dir(fullfile(fold,'*.png'));dir(fullfile(fold,'*.tif'))];

    tiffs = {listing(:).name};
    folds = {listing(:).folder};
    files = cellfun(@(x,y) fullfile(x,y),folds,tiffs,'UniformOutput',false);
    
    nameMol = zeros(1,length(listing));
    idxMov = zeros(1,length(listing));

    for i=1:length(listing) 
        spltName = strsplit(listing(i).name,'_mol_');
        nameMov = spltName{1};
        spltName2 = strsplit(spltName{2},'.');
        nameMol(i) = str2num(spltName2{1});
        names = cellfun(@(x) x.name,output,'un',false);
        idxMov(i) = find(cellfun(@(x) ~isempty(strfind(x,nameMov)),names));        
    end

% %            listing     
                
%     i=1;
    selected = [];
    for i=1:numImagesToShow:length(listing)-numImagesToShow+1
        try
            selected = [selected; (i-1)+select_image(listing,output,i,nameMol,idxMov,numImages(1),numImages(2))]; 
        end
    end
    if ~isempty(1:numImagesToShow:length(listing)-numImagesToShow+1)
        iLast = i+numImagesToShow;
    else
        iLast = 1;
    end
    
%     if isempty(i)
%         iLast = 1;
%     end
    if (iLast < length(files))
        try
            selected = [selected; iLast-1+select_image(listing,output,iLast,nameMol,idxMov,numImages(1),numImages(2))]; 
        end
    end
    
    allKymos = zeros(1,length(files));
    allKymos(selected) = 1;


    mkdir(foldOut,'good');
    mkdir(foldOut,'bad');
    % delete(fullfile(foldOut,'good/','*.tif'));
    % delete(fullfile(foldOut,'bad/','*.tif'));
    % 
    outputNew = cell(1,length(output));
    
    for i = 1:length(files)
        if allKymos(i) == 1
            copyfile(files{i},fullfile(foldOut,'good/',tiffs{i}))
        else
            copyfile(files{i},fullfile(foldOut,'bad/',tiffs{i}))
        end
    end
    
    goodBars = find(allKymos);
    
    outputNew = cell(1,length(outputNew));
    
    for i=1:length(output);
        curMols = find(idxMov(goodBars) == i);
        outputNew{i}.lineParams = output{i}.lineParams(curMols);
        outputNew{i}.xy = output{i}.xy(curMols);
        outputNew{i}.expBars = output{i}.expBars(curMols);
        %  outputNew{i}.expBars = delid;
        outputNew{i}.nanid = output{i}.nanid(curMols);
        outputNew{i}.idx = output{i}.idx(curMols);
        outputNew{i}.dotBars = output{i}.dotBars(curMols);
        outputNew{i}.boundaries = output{i}.boundaries(curMols);
        outputNew{i}.dots = output{i}.dots(curMols);
        outputNew{i}.trueedge = output{i}.trueedge(curMols);
        outputNew{i}.pos = output{i}.pos(curMols);
        outputNew{i}.name = output{i}.name;
        outputNew{i}.lengthLims = output{i}.lengthLims;
        outputNew{i}.widthLims = output{i}.widthLims;
        outputNew{i}.molScoreLim = output{i}.molScoreLim;
        outputNew{i}.dotScoreLim = output{i}.dotScoreLim;
        outputNew{i}.settings = output{i}.settings;
    end
    

 


end

