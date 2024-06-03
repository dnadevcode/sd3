function [images,names] = import_images_simple(experiment)
    % Import images from an experiment
    %
    %   Args:
    %       experiment - name of the file/folder
    %   Returns:
    %       images - loaded to mat
    %       names - names of images
tic;
folder = experiment.targetFolder;
fprintf('Analyzing data from the %s folder,\n',folder)
barFlag = experiment.barFlag;
dotFlag = experiment.dotFlag;
dotFlag2 = experiment.dotFlag2;

list = dir(folder);

if experiment.askForNumChannels
    numberOfChannels = str2num(questdlg('How many channels?', ...
	'How many channels', ...
	'1','2','3','2'));

else
    numberOfChannels = 2; % default number
end

if numberOfChannels==1
    barFlag = '';
end

if numberOfChannels~=3
    dotFlag2 = '';
end

images = {};
names = {};
for i = 1:numel(list)
  % Check if the file ends in 'C=1.tif'
    filename = fullfile(list(i).folder, list(i).name);
    [~,filen,ending] =fileparts(filename);
    switch ending(2:end)
        case 'czi' % todo: reduce try-catch
            % use bfopen to load czi
            try
                T = evalc(['data = bfopen(''', filename, ''');']);
            catch
                bfmatlabFold = uigetdir('pwd','Select folder with bfmatlab. ');
                addpath(genpath(bfmatlabFold));
                try
                    T = evalc(['data = bfopen(''', filename, ''');']);
                catch
                    warning('Failed to import czi file. Please check for bfmatlab');
                end
            end

            try % .czi only single frame
                images{end+1}.registeredIm = {double(data{1,1}{1,1})}; % Here just single frame
                [f1,f2,fend] = fileparts(filename);

                names{end+1} = f2;

                if ~isempty(barFlag) && numberOfChannels>=2
                    try
                        images{end}.dotIm = double(data{1,1}{2,1});
                    catch
                    end

                    if numberOfChannels>=3
                        try % if three channel
                            images{end}.dotIm2 = double(data{1,1}{3,1});
                        catch
                        end
                    end
                end
                data = [];
            catch
            end
        case 'tiff' % multitiff, only single frames
            images{end+1}.registeredIm{1} =   imread(filename,1);
            if ~isempty(barFlag)&& numberOfChannels>=2
                images{end}.dotIm =   imread(filename,2);
                 if numberOfChannels>=3
                    images{end}.dotIm2 =   imread(filename,3);
                 end
            end
            names{end+1} = filen;
        case 'tif'
            % loading a tiff
            isFile = isfile(filename);
            isBar = not(isempty(strfind(list(i).name, barFlag))) || isempty(barFlag);
            if isFile && isBar
                try
                    numMols = length(imfinfo(filename));
                    registeredIm = cell(1,numMols);
                    for idxTf=1:numMols
                        registeredIm{idxTf} =  imread(filename,idxTf);
                    end
                    images{end+1}.registeredIm = registeredIm;
                catch
                    fprintf('Could not interpret %s as an image.\n', list(i).name);
                    continue
                end
                
                names{end+1} = strrep(list(i).name, barFlag, ''); %sort away flag
                if not(isempty(barFlag))
                    dotPath = fullfile(list(i).folder, strrep(list(i).name, barFlag, dotFlag));
                    try
                        images{end}.dotIm = importdata(dotPath);
                    catch
                        fprintf('Did not find file %s.\n',dotPath);
                    end

                elseif isempty(barFlag) && isempty(dotFlag)
                    images{end}.dotIm = images{end}.registeredIm;
                end

                if not(isempty(dotFlag2))
                    dotPath2 = fullfile(list(i).folder, strrep(list(i).name, barFlag, dotFlag2));
                    try
                        images{end}.dotIm2 = importdata(dotPath2);
                    catch
                        fprintf('Did not find file %s.\n',dotPath2);
                    end
                end
            end
        
        otherwise
    end
end

if isempty(images)
%   throw(MException('image:import', 'No viable images found in target folder.'))
  warning(compose("No viable images found in target folder: %s", folder))
end


t=toc;
fprintf('Image import completed in %.1f seconds. \n',t);
