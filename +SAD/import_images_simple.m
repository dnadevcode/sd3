function [images,names] = import_images_simple(experiment,actions)
tic;
folder = experiment.targetFolder;
fprintf('Analyzing data from the %s folder,\n',folder)
barFlag = experiment.barFlag;
dotFlag = experiment.dotFlag;

list = dir(folder);
images = cell(1,floor(numel(list)/2));
imageSets = 0;
names = cell(1,floor(numel(list)/2));
for i = 1:numel(list)
    % Check if the file ends in 'C=1.tif'
    match = strfind(list(i).name,barFlag);
    if ~isempty(match)
        imageSets = imageSets + 1;
		names{imageSets} = list(i).name([1:match-1,match+length(barFlag):end]); %sort away flag
        images{imageSets}.registeredIm = cell(1);
		images{imageSets}.registeredIm{1} = importdata([folder,'/',list(i).name]);
        dotPath = [folder,'/',list(i).name(1:match-1),dotFlag,list(i).name(match+length(barFlag):end)];
            try
                images{imageSets}.dotIm = importdata(dotPath);
            catch
                fprintf('Did not find file %s.\n',dotPath);
            end
    end
    images(imageSets+1:end) = [];
	names(imageSets+1:end) = [];
end

%import PD.Core.Extraction.get_load_tiffs;


%[cutoutM] = get_load_tiffs(folder,actions.checkSameOrientation,actions.removeRegions,actions.makeSameSize);
 
%import PD.Core.Extraction.register_images_fast;
%images = register_images_fast(cutoutM);
t=toc;
fprintf('Image import completed in %.1f seconds. \n',t);
