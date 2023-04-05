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

list = dir(folder);

images = {};
names = {};
for i = 1:numel(list)
  % Check if the file ends in 'C=1.tif'
  filename = fullfile(list(i).folder, list(i).name);
  if isequal(filename(end-2:end),'czi')
      % use bfopen to load czi
      try
          T = evalc(['data = bfopen(''', filename, ''');']);
      catch
          bfmatlabFold = uigetdir('pwd','Select folder with bfmatlab');
          addpath(genpath(bfmatlabFold));
          try
              T = evalc(['data = bfopen(''', filename, ''');']);
          catch
              warning('Failed to import czi file');
          end

      end
        try
            images{end+1}.registeredIm = {double(data{1,1}{1,1})}; % Here just single frame
            [f1,f2,fend] = fileparts(filename);

            names{end+1} = f2; 
            try
              images{end}.dotIm = double(data{1,1}{2,1});
            catch
            end
            
            try % if three channel
              images{end}.dotIm2 = double(data{1,1}{3,1});
            catch
            end
            data = [];
        catch
        end

  else
      % loading a tiff
  isFile = isfile(filename);
  isBar = not(isempty(strfind(list(i).name, barFlag))) || isempty(barFlag);
  if isFile && isBar
    try
        numMols = length(imfinfo(filename));
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
  end
  end
end

if isempty(images)
%   throw(MException('image:import', 'No viable images found in target folder.'))
  warning(compose("No viable images found in target folder: %s", folder))
end

%import PD.Core.Extraction.get_load_tiffs;


%[cutoutM] = get_load_tiffs(folder,actions.checkSameOrientation,actions.removeRegions,actions.makeSameSize);

%import PD.Core.Extraction.register_images_fast;
%images = register_images_fast(cutoutM);
t=toc;
fprintf('Image import completed in %.1f seconds. \n',t);
