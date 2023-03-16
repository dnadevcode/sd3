function [barcodes, dotScoreMin] = sdd_extract_barcodes(movies, sets, lengthLims, imageNumber, runNo, bgPixels,tiles)
    % sdd_extract_barcodes
    %
    %   Args:
    %         movies - movie structure
    %
    %   Returns:
    %       barcodes - detected barcodes
    %       dotScoreMin - minimum dot score threshold

%%%%%%%%%%%%%%%%% TWEAK PARAMETER %%%%%%%%%%%%%%%%%
  sPer = round(sets.sigma); % Number of pixel width of molecule "kymograph"
  %%%%%%%%%%%%%%%%% TWEAK PARAMETER %%%%%%%%%%%%%%%%%

  extractionMethod = sets.extractionMethod;
  % Convert images to kymographs
  [kymos, barcodes.lineParams, barcodes.xy,barcodes.distance] = get_kymos_from_movies(movies.molM, movies.bwM, sPer, extractionMethod);

%%  
%     figure,
%     idx = 1;
%     mol = movies.molM{idx};
%     xF = barcodes.xy{idx}{1};
%     yF = barcodes.xy{idx}{2};
%     
%         tiledlayout(1,2);nexttile
%         imagesc(mol')
%           nexttile
%         imagesc(mol')
%             hold on
% %             
%             plot(xF,yF,'redx')
% % %             
%         colormap(gray)
      %%  
        
% %   % plot
%     idx=1
%     figure,
%     tiledlayout(1,2);nexttile
%     if extractionMethod == 2
% %         imagesc(movies.molM{idx}');
%     else
%         imagesc(movies.molM{idx})
%         nexttile
%         imagesc(movies.dotM{idx})
%         hold on
%         plot(-barcodes.lineParams{idx}(1)*(1:size(movies.molM{idx},2))+barcodes.lineParams{idx}(2),'redx')
%     end
%     colormap(gray)

%     %   
%     imwrite(uint16(movies.dotM{idx}),'ex.tif')
%     sPer = 0;
    
  % Extract barcodes from kymographs
    [barcodes.expBars, barcodes.expStats, barcodes.delid, barcodes.nanid] = trans_bar_extr(kymos, sets);
    barcodes.idx = 1:length(movies.molM);
    barcodes.idx( barcodes.delid) = [];
  if isfield(movies, 'dotM')
      

      if   extractionMethod ==1
    % Extract dotBarcodes from their images
        [barcodes.dotBars,barcodes.boundaries] = cellfun(@(x,y) get_dot_kymo(x, y(1) , y(2), sPer) ,movies.dotM, barcodes.lineParams,'un',false);
        barcodes.dotBars(barcodes.delid) = [];

      else
    % curved molecule
        barcodes.xy(barcodes.delid) = []; % so it doesn't need to compute those
        movies.dotM(barcodes.delid) = [];
        [barcodes.dotBars,barcodes.boundaries] = cellfun(@(x,y) get_curved_kymo(x, y{1} , y{2}, sPer), movies.dotM, barcodes.xy,'un',false);
      end
      

%     [barcodes.dotBars, lineParams2] = get_kymos_from_movies(movies.dotM, movies.bwM, sPer);
%     barcodes.xy( barcodes.delid) = [];
%     barcodes.boundaries( barcodes.delid) = [];
%     barcodes.lineParams( barcodes.delid) = [];

    % todo: SFW detection here - but for this need to remove noise first
    
%     if extractionMethod == 2
%         numavgPts = length(round(-669/509 * sets.sigma):round(669/509 * sets.sigma));
%         bgPixels = movmean(bgPixels,  numavgPts); % but for method 2, we take a single point?
%     end

    % old:
    [barcodes.dots, dotScoreMin] = sdd_detect_dots(barcodes.dotBars, sets, imageNumber, movies.imageName, bgPixels,tiles);
%     if sets.dotDet2D % detect 2d dots      
%         % for each point need to find 
%         %%
% %         [barcodes.dots, dotScoreMin] = sdd_detect_dots_2D(movies.dotM, optics, sets.dotMargin, sets.dotScoreMin, lengthLims, imageNumber, movies.imageName, sets, bgPixels,tiles);
%     end
  else
    dotScoreMin = 'NA';
  end
%%  also for spline
%     idx = 2;
%     figure
%     imagesc(movies.dotM{idx})
% %         imagesc(movies.molM{idx})
% % 
%     hold on
%     plot(barcodes.boundaries{idx}(3):barcodes.boundaries{idx}(4),-barcodes.lineParams{idx}(1)*(barcodes.boundaries{idx}(3):barcodes.boundaries{idx}(4))+barcodes.lineParams{idx}(2),'red')
%     angle = atan(barcodes.lineParams{idx}(1));
% % 
%     vOff = barcodes.boundaries{idx}(1);
%     hOff = barcodes.boundaries{idx}(3);
% 
%     for j = 1:numel(barcodes.dots{idx}.locations)
%             dy = -sin(angle)*(barcodes.dots{idx}.locations(j)-1+barcodes.dots{idx}.leftOffset);%-sin(angle)*barcodes.dots{idx}.leftOffset); % this ofset should be from (1,1) of the map
%             dx = cos(angle)*(barcodes.dots{idx}.locations(j)-1+barcodes.dots{idx}.leftOffset);
%             y = vOff+dy-1;%barcodes.lineParams{molIdx}(2)+dy-1;
%             x = hOff+dx;
%             plot(x,y,'gx'); % maybe too much to plot also this info
%             %         str = sprintf('I = %.1f, depth = %.1f',barcodes.dots{molIdx}.val(j),barcodes.dots{molIdx}.depth(j));
%             %         text(x-5,y-5,str,'Color','white');
%     end
%     colormap(gray)
%%
  if sets.saveBars
    addpath(genpath( subsref(dir(sets.targetFolder), substruct('.', 'folder'))));

    saveIdx = 0;

    for i = 1:length(barcodes.expBars)
      saveIdx = saveIdx + 1;

      while sum(barcodes.delid == saveIdx) == 1
        saveIdx = saveIdx + 1;
      end

      bc = barcodes.expBars{i}.rawBarcode;
      saveName = regexprep(movies.imageName, '[.]\S{1,4}', '');
      folderName = subsref(dir(sets.targetFolder), substruct('.', 'folder'));
      imwrite(uint16(bc), fullfile(folderName, ['barcodes_run', num2str(runNo)], [saveName, 'barcode', num2str(saveIdx), '.tif']));
      
      
      km = barcodes.expBars{i}.kymo;
      imwrite(uint16(km), fullfile(folderName, ['kymos_run', num2str(runNo)], [saveName, 'kymo', num2str(saveIdx), '.tif']));

      if isfield(movies, 'dotM')
        db = barcodes.dotBars{i};
        db = db(~isnan(db));
        imwrite(uint16(db), fullfile(folderName, ['dotbars_run', num2str(runNo)], [saveName, 'dotbar', num2str(saveIdx), '.tif']));
      end

    end

  end
