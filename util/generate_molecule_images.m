function [molM,bwM,pos] = generate_molecule_images(D,registeredIm,saveName,edgePx,actions)
cellLength = length(D);
pos = cell(1,cellLength);
molM = cell(1,cellLength);
bwM = cell(1,cellLength);
mkdir([saveName,'/molsandbars']);

for k = 1:cellLength
    yMin = max(1,min(D{k}(:,1))-edgePx);
    yMax = min(edgePx+max(D{k}(:,1)),size(registeredIm{1},1));
    xMin = max(1,min(D{k}(:,2))-edgePx);
    xMax = min(edgePx+max(D{k}(:,2)),size(registeredIm{1},2));
    molMov = nan(yMax-yMin+1,xMax-xMin+1,length(registeredIm));
    bwMov = nan(yMax-yMin+1,xMax-xMin+1);
	 for i = yMin:yMax
	 	  for j = xMin:xMax
				in = inpolygon(j,i,D{k}(:,2),D{k}(:,1));
				if in
					for l = 1:length(registeredIm)
						molMov(i-yMin+1,j-xMin+1,l) = registeredIm{l}(i,j);
					end
					bwMov(i-yMin+1,j-xMin+1) = 1; 
				end
		  end
	 end
	 if actions.saveMolecules
		 molSaveName = [saveName,'/molsandbars/mol',num2str(k),'mov.tiff'];
		 imwrite(uint16(molMov(:,:,1)),molSaveName);
		 for i = 2:length(registeredIm)
   	     imwrite(uint16(molMov(:,:,i)),molSaveName,'writemode','append');
		 end
%		 bwSaveName = [saveName,'/bw',num2str(k),'.tiff'];	 
%		 imwrite(uint16(bwMov(:,:,1)),bwSaveName);
	 end
	 pos{k} = D{k}(1,:);
  	 molM{k} = molMov;
	 bwM{k} = bwMov;
end

molM = molM(~cellfun('isempty',molM));
bwM = bwM(~cellfun('isempty',bwM));
