function [newBars,barStats,delId,nanid] = trans_bar_extr(kymos,optics,actions,sets)

edgeLen = round(sets.deltaCut * optics.sigma);

% A more transparent alignment/barcode extraction routine
nmol = length(kymos);
newBars = cell(1,nmol);
nanid = zeros(1,nmol);
for i = 1:nmol
  initbar = nanmean(kymos{i},1);
  newBars{i}.rawBarcode = initbar(~isnan(initbar));
  nanid(i) = find(~isnan(initbar),1,'first');
  bitmaskLen = length(newBars{i}.rawBarcode);
  newBars{i}.rawBitmask = ones(1,bitmaskLen);
  newBars{i}.rawBitmask([1:min(edgeLen,bitmaskLen),max(bitmaskLen - edgeLen + 1,1):end]) = 0;
end

% if actions.getConsensus == 1
% 	commonLength = ceil(mean(cellfun(@(x) length(x.rawBarcode),newBars)));
% 	barStats.commonLength = commonLength;
% 	for i = 1:nmol
% 		molLen = length(newBars{i}.rawBarcode);
% 		if molLen > 1
% 			v = linspace(1,molLen,commonLength);
% 			newBars{i}.stretchedBarcode = interp1(newBars{i}.rawBarcode,v);
% 			newBars{i}.stretchedBitmask = newBars{i}.rawBitmask(round(v));
% 		end
% 	end
% end

% Remove barcodes with less than 3 trusted, numerical points
delId = [];
for i = 1:nmol
  detail = sum(~isnan(newBars{i}.rawBarcode).*newBars{i}.rawBitmask);
  if detail < 2*edgeLen+3 % Ensure that there is always 3 data points to compare between theory and exp
    delId = [delId i];
  end
end
newBars(delId) = [];
nanid(delId) = [];
fprintf('%i viable barcodes found.\n',length(newBars));

barStats.lengthAverage = mean(cellfun(@(x) length(x.rawBarcode),newBars));
barStats.lengthStd = std(cellfun(@(x) length(x.rawBarcode),newBars));
