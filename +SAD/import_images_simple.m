function [images,names] = import_images_simple(experiment)
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
  isFile = isfile(filename);
  isBar = not(isempty(strfind(list(i).name, barFlag))) || isempty(barFlag);
  if isFile && isBar
    try
      images{end+1}.registeredIm = imread(filename);
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

if isempty(images)
  throw(MException('image:import', 'No viable images found in target folder.'))
end

%import PD.Core.Extraction.get_load_tiffs;


%[cutoutM] = get_load_tiffs(folder,actions.checkSameOrientation,actions.removeRegions,actions.makeSameSize);

%import PD.Core.Extraction.register_images_fast;
%images = register_images_fast(cutoutM);
t=toc;
fprintf('Image import completed in %.1f seconds. \n',t);
