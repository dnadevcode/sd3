function dnarec_folder_scan(target, sets)

if nargin == 0
  target = '/testfolder/';
end

% Recursion:
% Make list of possible folders
% check if any of them have subfolders or are output
% if not - add to list of data folders and remove from list of possibles
% if output - remove from list of possibles
% if subfolders exist - add each to list of possible - remove mother
% folder
% terminate when list is emptys

dataFolders = search_folder(target);
if nargin < 2
  for i = 1:numel(dataFolders)
    dnarec_skel(dataFolders{i});
  end
else
  for i = 1:numel(dataFolders)
    dnarec_skel(dataFolders{i}, sets);
  end
end

end

function dataFolders = search_folder(target)
possibles{1} = target;
nDataFolders = 0;
while numel(possibles) > 0
  for i = 1:numel(possibles)
    isOutput = filter_output(possibles{i});
    [hasSubfolder,subfolders] = filter_master(possibles{i});
    isLink = filter_link(possibles{i});
    isEmpty = not(hasSubfolder) && not(filter_images(possibles{i}));
    %Check if possibles is an output folder
    if isOutput || isLink || isEmpty
      % Remove element from possibles
    else
      if hasSubfolder
        % Add each subfolder to possibles
        possibles = [possibles subfolders];
      end
      nDataFolders = nDataFolders + 1;
      % Add element to datafolders
      dataFolders{nDataFolders} = possibles{i};
    end
    possibles{i} = [];
  end
  % Remove empty entries from possibles
  emptyEntries = cellfun(@(x) isempty(x),possibles);
  possibles(emptyEntries) = [];
end
end

function isOutput = filter_output(folder)
isOutput = contains(folder,'molsandbars') | contains(folder,'movies');
end

function hasImages = filter_images(folder)
content = dir(folder);
fContent = content(not([content.isdir]));
hasImages = false;
for f=fContent'
  try
    imfinfo(fullfile(f.folder, f.name));
    hasImages = true;
    return
  catch
    continue
  end
end
end

function [hasSubfolder,subfolders] = filter_master(folder)
content = dir(folder);
subfolders = cell(1,numel(content));
nSubfolders = 0;
isFolder = [content.isdir];
content = content(isFolder);
for i = 1:numel(content)
  isOutput = filter_output(content(i).name);
  isLink = filter_link(content(i).name);
  if ~(isLink || isOutput)
    nSubfolders = nSubfolders + 1;
    subfolders{nSubfolders} = fullfile(content(i).folder, content(i).name);
  end
end
hasSubfolder = nSubfolders > 0;
end

function isLink = filter_link(folder)
isLink = strcmp(folder,'.') | strcmp(folder,'..');
end
