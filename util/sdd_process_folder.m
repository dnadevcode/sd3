function [output,hPanelResult,images,movies,barcodes] = sdd_process_folder(dataFold, sets, tsHCC)
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

    if nargin < 3
        % use no gui.
       tsHCC = [];
    end

    if nargout >=3
        images = [];
        movies = [];
        barcodes = [];
    end
    
    
    hPanelResult = [];



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
    kymoName = fullfile(folderName, ['kymos_run', num2str(runNo)]);
    dotName = fullfile(folderName, ['dotbars_run', num2str(runNo)]);
    molName = fullfile(folderName, ['molecules_run', num2str(runNo)]);
    resName = fullfile(folderName, ['results_run', num2str(runNo), '.txt']);
    barcodeExist = isfolder(barName);
    dotbarExist = isfolder(dotName);
    moleculesExist = isfolder(molName);
    resExist = isfile(resName);
    kymoExist = isfolder(kymoName);
    outputExist = barcodeExist || dotbarExist || moleculesExist || resExist || kymoExist;
  end

  if sets.saveMolecules
    mkdir(molName)
  end

  if sets.saveBars
    mkdir(kymoName)
    mkdir(barName)
    mkdir(dotName)
  end
  
    import SAD.denoise_images;
    hPanelResult = cell(1,numel(images));
    btn = cell(1,numel(images)); % helping plot button
    btn2 = cell(1,numel(images)); % helping plot button

  for i = 1:numel(images)
    fprintf('\nAnalysing image %s.\n', imageNames{i});
    if ~isempty(tsHCC)

    hPanelResult{i} = uitab(tsHCC, 'title',strcat([ imageNames{i} ' results']));
    
    if sets.showMolecules || sets.showScores
        % create tabs to display results from the analysis
        hTabgroup = uitabgroup('Parent',hPanelResult{i});

        hResScores= uitab(hTabgroup, 'title',strcat('Scores'));
        t = tiledlayout(hResScores,2,2,'TileSpacing','compact','Padding','compact');
        tiles.molScores = nexttile(t);
        tiles.dotScores = nexttile(t);
        tiles.dotScoresFilt = nexttile(t);
        tiles.bgScores = nexttile(t);

        hResplot = uitab(hTabgroup, 'title',strcat('Detected molecules'));
        t = tiledlayout(hResplot,1,2,'TileSpacing','tight','Padding','tight');
        tiles.molDet = nexttile(t);
        tiles.dotDet = nexttile(t);
        linkaxes([ tiles.molDet tiles.dotDet  ])

        hResFilt = uitab(hTabgroup, 'title',strcat('Filtered'));
        t = tiledlayout(hResFilt,1,2,'TileSpacing','tight','Padding','tight');
        tiles.logFilt = nexttile(t);
        tiles.bg = nexttile(t);
    else
        tiles = [];
    end
    else
        tiles = [];
  end
    % Detect centers (no denoising)
    cleanImages.registeredIm = images{i}.registeredIm;
    cleanImages.imAverage =  averageImages(images{i}.registeredIm);
    cleanImages.centers = [];
    sets.denoiseImages = 0;% not used at the moment
    if sets.denoiseImages
        cleanImages = denoise_images(images{i}.registeredIm,sets);
    end
    
    if isfield(images{i}, 'dotIm')
        cleanImages.dotIm = images{i}.dotIm;
        if sets.denoiseDotImages
            SE = strel("rectangle",[20 20]);
            cleanImages.dotIm  = imtophat(double(cleanImages.dotIm ),SE);
        end
    end

    % Segment image
    [movies, ~, sets, lengthLims, molScoreLim, widthLims, bgPixels, bgPixels2] = sdd_segment_image(cleanImages, imageNames{i}, i, runNo,sets,tiles);

    if isempty(bgPixels2)
        bgPixels2 = bgPixels;
    end
    
%     sets.extractionMethod = 2;
    % Extract barcodes
    [barcodes, dotScoreLim] = sdd_extract_barcodes(movies, sets, lengthLims, i, runNo, bgPixels2, tiles);

%   % plot molecule with mask
%       if sets.saveMolecules
% % 
%           idx=find(barcodes.idx==5)
%           figure,imagesc(movies.dotM{barcodes.idx(idx)})
%     %         imagesc(movies.molM{barcodes.idx(idx)})
% 
%           hold on
%           plot(barcodes.xy{barcodes.idx(idx)}{2},barcodes.xy{barcodes.idx(idx)}{1},'redx')
% %           pos = barcodes.dots{idx}.locations+barcodes.nanid(idx);
% %           plot(barcodes.xy{idx}{2}(pos),barcodes.xy{idx}{1}(pos),'blackx')
%     %             plot(-barcodes.lineParams{idx}(1)*(1:size(movies.molM{idx},2))+barcodes.lineParams{idx}(2),'redx')
% 
% %           plot(-barcodes.lineParams{idx}(1)*(1:size(movies.molM{idx},2))+barcodes.lineParams{idx}(2),'redx')
%     % %     %   
%         colormap(gray)
% % % %     imwrite(uint16(movies.dotM{idx}),'ex.tif')
% % % %     sPer = 0;
%       end

    % Calculate p-values for specific sequence
    output{i} = barcodes;
    if sets.showMolecules % optional
        output{i}.movies = movies; 
    end
    output{i}.trueedge = movies.trueedge;
    output{i}.pos = movies.pos;

    output{i}.name = imageNames{i};
    output{i}.lengthLims = lengthLims;
    output{i}.widthLims = widthLims;
    output{i}.molScoreLim = molScoreLim;
    output{i}.dotScoreLim = dotScoreLim;
    output{i}.median = median(cleanImages.imAverage(:));
    output{i}.settings = sets;
    output{i}.molRunFold = movies.molRunName;
    output{i}.runNo = movies.runNo;
    
    % csv print:
  import SAD.csv_print
  resultsName = csv_print(output, sets, runNo, i);

      
    
    if sets.showMolecules
      % Mark barcodes in molecule image
        [output{i}.dotLocsGlobal] = sdd_mark_bars(  output{i}.movies, output{i},tiles,sets.extractionMethod);
        if i==numel(images) % define only for last image
            tb = axtoolbar(tiles.dotDet , 'default');
            btn{i} = axtoolbarbtn(tb,'Icon',1+63*(eye(25)),'Tooltip','Detailed molecule plot');
            btn{i}.ButtonPushedFcn = @callbackDetailedPlot; 
            btn2{i} = axtoolbarbtn(tb,'Icon',1+63*(triu(ones(25))),'Tooltip','Zoomed-in molecule plot 2');
            btn2{i}.ButtonPushedFcn = @callbackDetailedPlot2; 
        end
    else
        [output{i}.dotLocsGlobal] = dot_locs_global( movies, barcodes ,sets.extractionMethod);
    end
    
  end

%   fig3(images, movies,barcodes)

  import SAD.dnarec_print
  resultsName = dnarec_print(output, sets, runNo);
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

  
  
    
    function callbackDetailedPlot(src, event)
        answer = inputdlg('Select molecule(s) to analyse in detail:',...
             'Mol analysis', [1 50]);
        user_val = str2num(answer{1});
        import Core.AnalysisPlot.detailed_analysis_plot;
        detailed_analysis_plot(movies,barcodes,user_val)

    end

    function callbackDetailedPlot2(src, event)
        answer = inputdlg('Select molecule(s) to analyse in detail:',...
             'Mol analysis', [1 50]);
        user_val = str2num(answer{1});
        import Core.AnalysisPlot.detailed_analysis_plot_nicefig;
        detailed_analysis_plot_nicefig(images, movies,barcodes,user_val)

    end
end
