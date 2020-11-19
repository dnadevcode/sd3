function [barcodes,dotScoreMin] = extract_barcodes(movies,optics,lengthLims,imageNumber,runNo,sets,experiment,actions)

%import PD.Core.Extraction.get_kymos_from_movies;
%%%%%%%%%%%%%%%%% TWEAK PARAMETER %%%%%%%%%%%%%%%%%
sPer = round(optics.sigma); % Number of pixel width of molecule "kymograph"
%%%%%%%%%%%%%%%%% TWEAK PARAMETER %%%%%%%%%%%%%%%%%

% Convert images to kymographs
[kymos,barcodes.lineParams] = get_kymos_from_movies(movies.molM,movies.bwM,sPer);

% Extract barcodes from kymographs
[barcodes.expBars,barcodes.expStats,barcodes.delid] = trans_bar_extr(kymos,optics,actions,sets);

if isfield(movies,'dotM')
  % Extract dotBarcodes from their images
  [barcodes.dotBars,~] = get_kymos_from_movies(movies.dotM,movies.bwM,sPer);
  barcodes.dotBars(barcodes.delid) = [];
  [barcodes.dots,dotScoreMin] = detect_dots(barcodes.dotBars,optics,sets.dotMargin,experiment.dotScoreMin,lengthLims,imageNumber,movies.imageName,actions);
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
    saveName = regexprep(movies.imageName, '[.]\S{1,4}', '');
    folderName = subsref(dir(experiment.targetFolder), substruct('.', 'folder'));
    imwrite(uint16(bc), fullfile(folderName, ['barcodes_run',num2str(runNo)], [saveName,'barcode',num2str(saveIdx),'.tif']));
    if isfield(movies,'dotM')
      db = barcodes.dotBars{i};
      db = db(~isnan(db));
      imwrite(uint16(db), fullfile(folderName, ['dotbars_run',num2str(runNo)], [saveName,'dotbar',num2str(saveIdx),'.tif']));
    end
  end
end
