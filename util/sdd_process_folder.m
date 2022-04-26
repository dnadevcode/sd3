function [output,hPanelResult] = sdd_process_folder(dataFold, sets, tsHCC)
    % sdd_process_folder
    % This routine analyses the image contents of the folder specified in the "exp" struct
    % via the functions specified in the "functions" struct according to te actions
    % specified in the "actions" struct.
    %
    %   Args:
    %       dataFold, sets,tsHCC
    %
    %   Returns:
    %       output, hPanelResult
    %
    
    
    hPanelResult = [];


  mFilePath = mfilename('fullpath');
  mfolders = split(mFilePath, {'\', '/'});
  utilPath = fullfile(mfolders{1:end - 2}, 'util');

  if strcmp(mFilePath(1), '/')
    utilPath = strcat('/', utilPath);
  end

  addpath(utilPath);
  [~, lwid] = lastwarn;

  if strcmp(lwid, 'MATLAB:mpath:nameNonexistentOrNotADirectory')
    error('Unexpected error when asserting source folder path.')
  end

  % addpath('util/lldev/src/MATLAB/');
  import SAD.dnarec_print

  % Perform preliminary check to see if all required files are accessible
  % check = prel_check(experiment,functions,actions);
  %
  % if check
  % Import and crop images
    sets.targetFolder = dataFold;
    import SAD.import_images_simple;
    [images, imageNames] = import_images_simple(sets);
%   [images, imageNames] = functions.img_import();
  output = cell(1, numel(images));

  if isempty(images) || prompt_figure_excess(length(images), sum(cellfun(@(x) isfield(x, 'dotIm'), images)), sets)
    return
  end

  % Specify save name
  % slashes = strfind(experiment.targetFolder,'/');
  % folderName = experiment.targetFolder(slashes(end-1)+1:slashes(end)-1);
  folderName = subsref(dir(sets.targetFolder), substruct('.', 'folder'));
  runNo = 0;
  outputExist = 1;

  while outputExist
    runNo = runNo + 1;
    %     if runNo > 98
    %         fprintf('Only %i run outputs supported, overwriting last result (run %i).\n',runNo,runNo);
    %         break;
    %     end
    barName = fullfile(folderName, ['barcodes_run', num2str(runNo)]);
    dotName = fullfile(folderName, ['dotbars_run', num2str(runNo)]);
    molName = fullfile(folderName, ['molecules_run', num2str(runNo)]);
    resName = fullfile(folderName, ['results_run', num2str(runNo), '.txt']);
    barcodeExist = isfolder(barName);
    dotbarExist = isfolder(dotName);
    moleculesExist = isfolder(molName);
    resExist = isfile(resName);
    outputExist = barcodeExist || dotbarExist || moleculesExist || resExist;
  end
  


  if sets.saveMolecules
    mkdir(molName)
  end

  if sets.saveBars
    mkdir(barName)
    mkdir(dotName)
  end
  
    import SAD.denoise_images;
    hPanelResult = cell(1,numel(images));

  for i = 1:numel(images)
    fprintf('\nAnalysing image %s.\n', imageNames{i});

    hPanelResult{i} = uitab(tsHCC, 'title',strcat([ imageNames{i} ' results']));
    
    if sets.showMolecules || sets.showScores
        % create tabs to display results from the analysis
        hTabgroup = uitabgroup('Parent',hPanelResult{i});

        hResScores= uitab(hTabgroup, 'title',strcat('Scores'));
        t = tiledlayout(hResScores,1,2,'TileSpacing','tight','Padding','tight');
        tiles.molScores = nexttile(t);
        tiles.dotScores = nexttile(t);
        hResplot = uitab(hTabgroup, 'title',strcat('Detected molecules'));
        t = tiledlayout(hResplot,1,2,'TileSpacing','tight','Padding','tight');
        tiles.molDet = nexttile(t);
        tiles.dotDet = nexttile(t);
        linkaxes([ tiles.molDet tiles.dotDet  ])

        hResFilt = uitab(hTabgroup, 'title',strcat('Filtered'));
        t = tiledlayout(hResFilt,1,2,'TileSpacing','tight','Padding','tight');
        tiles.logFilt = nexttile(t);

%         t = tiledlayout(hPanelResult,4,2,'TileSpacing','tight','Padding','tight');
%         tiles.molScores = nexttile(t);
%         tiles.dotScores = nexttile(t);
%         tiles.molDet = nexttile(t,[2 1]);
%         tiles.dotDet = nexttile(t,[2 1]);
%         linkaxes([ tiles.molDet tiles.dotDet  ])
%         tiles.logFilt = nexttile(t);

%     if sets.showMolecules %create tile layout if plotting the results
% %         t = tiledlayout(hPanelResult,3,2,'TileSpacing','tight','Padding','tight');
%     else
    else
        tiles = [];
    end
 
    % Denoise images and remove artifacts
    cleanImages = denoise_images(images{i}.registeredIm,sets);

    if isfield(images{i}, 'dotIm')
      cleanImages.dotIm = images{i}.dotIm;
    end

    % Segment image
    [movies, ~, optics, lengthLims, molScoreLim, widthLims, bgPixels, bgPixels2] = sdd_segment_image(cleanImages, imageNames{i}, i, runNo,sets,tiles);

    if isempty(bgPixels2)
        bgPixels2 = bgPixels;
    end
    % Extract barcodes
    [barcodes, dotScoreLim] = sdd_extract_barcodes(movies, optics, lengthLims, i, runNo, bgPixels2, sets, tiles);

%   % plot
%       idx=16
%       figure,imagesc(movies.molM{idx})
%       hold on
%       plot(-barcodes.lineParams{idx}(1)*(1:size(movies.molM{idx},2))+barcodes.lineParams{idx}(2),'redx')
% %     %   
%     imwrite(uint16(movies.dotM{idx}),'ex.tif')
%     sPer = 0;
    

    if sets.showMolecules
      % Mark barcodes in molecule image
      sdd_mark_bars(movies, barcodes,tiles);
    end

    % Calculate p-values for specific sequence
    %	output = functions.bc_analyse(barcodes,optics);i
    output{i} = barcodes;
    output{i}.name = imageNames{i};
    output{i}.lengthLims = lengthLims;
    output{i}.widthLims = widthLims;
    output{i}.molScoreLim = molScoreLim;
    output{i}.dotScoreLim = dotScoreLim;
    output{i}.median = median(cleanImages.imAverage(:));
  end

  import SAD.dnarec_print
  resultsName = dnarec_print(output, sets, sets, optics, runNo, sets);
  fprintf('\n-------------------------------------------------------------------\n');
  fprintf('Analysis complete\n');
  fprintf('Results saved in %s', resultsName);
  fprintf('\n-------------------------------------------------------------------\n');

  % Run the consensus barcode routine on the barcodes
  %	consensusStruct = functions.consensus_gen(barcodes);

  % Print results
  %	dnarec_print(output)
  % Mark significants
  %	mark_matches(output,movies.pos,barcodes.delid,0,0.09)

  %	% Save results
  %	dnarec_save(images,cleanImages,movies,scores,barcodes,output,experiment,consensusStruct);
  % else
  % 	fprintf('Preliminary check failed. Analysis aborted.\n');
  % 	output = NaN;
end
