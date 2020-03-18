function [barcodes,dotScoreMin] = extract_barcodes(movies,optics,lengthLims,imageNumber,runNo,sets,experiment,actions)

%import PD.Core.Extraction.get_kymos_from_movies;
%%%%%%%%%%%%%%%%% TWEAK PARAMETER %%%%%%%%%%%%%%%%%
sPer = optics.sigma;
%%%%%%%%%%%%%%%%% TWEAK PARAMETER %%%%%%%%%%%%%%%%%

% Convert images to kymographs
[kymos,barcodes.lineParams] = get_kymos_from_movies(movies.molM,movies.bwM,sPer);

% Extract barcodes from kymographs
[barcodes.expBars,barcodes.expStats,barcodes.delid] = trans_bar_extr(kymos,optics,actions,sets);

if isfield(movies,'dotM')
	% Extract dotBarcodes from their images
	[barcodes.dotBars,~] = get_kymos_from_movies(movies.dotM,movies.bwM,sPer);
    barcodes.dotBars(barcodes.delid) = [];
	[barcodes.dots,dotScoreMin] = detect_dots(barcodes.dotBars,optics,sets.dotMargin,lengthLims,imageNumber,movies.imageName,actions);
else
    dotScoreMin = 'NA';
end

if actions.saveBars
	saveIdx = 0;
	for i = 1:length(barcodes.expBars)
		saveIdx = saveIdx + 1;
		while sum(barcodes.delid == saveIdx) == 1
			saveIdx = saveIdx + 1;
		end
		bc = barcodes.expBars{i}.rawBarcode;
		imwrite(uint16(bc),[experiment.targetFolder,'barcodes_run',num2str(runNo),'/',movies.imageName(1:end-4),'barcode',num2str(saveIdx),'.tif']);
		if isfield(movies,'dotM')
			db = barcodes.dotBars{i};
			db = db(~isnan(db));
			imwrite(uint16(db),[experiment.targetFolder,'dotbars_run',num2str(runNo),'/',movies.imageName(1:end-4),'dotbar',num2str(saveIdx),'.tif']);
		end
	end
end
