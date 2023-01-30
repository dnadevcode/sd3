function [molM,bwM,dotM,pos] = generate_molecule_images_fast(D,L,registeredIm,dotIm,folder,runNo,imageName,edgePx,actions)
cellLength = length(D);
pos = cell(1,cellLength);
molM = cell(1,cellLength);
bwM = cell(1,cellLength);
dotM = cell(1,cellLength);

% Fix to simplify surrounding code.
if not(iscell(registeredIm))
  registeredIm = {registeredIm};
end

for k = 1:cellLength
  yMin = max(1,min(D{k}(:,1))-edgePx);
  yMax = min(edgePx+max(D{k}(:,1)),size(registeredIm{1},1));
  xMin = max(1,min(D{k}(:,2))-edgePx);
  xMax = min(edgePx+max(D{k}(:,2)),size(registeredIm{1},2));
  dY = yMax-yMin+1;
  dX = xMax-xMin+1;
  molMov = nan(dY,dX,length(registeredIm));
  bwMov = nan(dY,dX);
  
  labelCut = L(yMin:yMax,xMin:xMax);
  mask = double(labelCut == k);
  for l = 1:length(registeredIm)
    regCut = registeredIm{l}(yMin:yMax,xMin:xMax);
    molMov(:,:,l) = double(regCut).*mask;
    molMov(molMov == 0) = NaN;
  end
  if isempty(dotIm)
    dotMov = [];
  else
    dotCut = dotIm(yMin:yMax,xMin:xMax);
    dotMov = double(dotCut).*mask;
    dotMov(dotMov == 0) = NaN;
  end
  bwMov(mask == 1) = 1;
  
  if actions.saveMolecules
    [~,de,~] = fileparts(imageName);
%     de= strcat(de,['_mov=_' num2str(k) '_']);
    molSaveName = fullfile(folder, ['molecules_run',num2str(runNo)], [de,'_mol_',num2str(k),'.tif']);

    imwrite(uint16(regCut(:,:,1)),molSaveName);

    if ~isempty(dotCut)
      imwrite(uint16(dotCut),molSaveName,'writemode','append');
    end
    
    if ~isempty(dotMov)
      imwrite(uint16(dotMov),molSaveName,'writemode','append');
    end
        imwrite(uint16(bwMov(:,:,1)),molSaveName,'writemode','append');
        imwrite(uint16(molMov(:,:,1)),molSaveName,'writemode','append');
    for i = 2:length(registeredIm)
        imwrite(uint16(molMov(:,:,i)),molSaveName,'writemode','append');
    end
  end
  %	 pos{k} = D{k}(1,:);
  pos{k} = [yMin xMin];
  molM{k} = molMov;
  bwM{k} = bwMov;
  dotM{k} = dotMov;
end

molM = molM(~cellfun('isempty',molM));
bwM = bwM(~cellfun('isempty',bwM));
if isempty(dotIm)
  dotM = [];
else
  dotM = dotM(~cellfun('isempty',dotM));
end
