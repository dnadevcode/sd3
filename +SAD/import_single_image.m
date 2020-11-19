function [images,names] = import_single_image(targetFile,experiment,actions)

images = cell(1);
names = cell(1);

names{1} = targetFile(1:end-7);
folder = experiment.targetFolder;
images{1}.registeredIm = cell(1);

impath = [folder,names{1},'C=1.tif'];
try
  images{1}.registeredIm{1} = importdata(impath);
catch
  fprintf('Did not find file %s.\n',impath);
end
dotpath = [folder,names{1},'C=0.tif'];
try
  images{1}.dotIm = importdata(dotpath);
catch
  fprintf('Did not find file %s.\n',dotpath);
end
