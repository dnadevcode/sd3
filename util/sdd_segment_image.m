function [movies, scores, optics, lengthLims, lowLim, widthLims, bgCutOut, bgCutOut2] = sdd_segment_image(images, imageName, imageNumber, runNo, sets, tiles)

  registeredIm = images.registeredIm;
  imAverage = images.imAverage;
  %imDenoised = images.imDenoised;
  centers = images.centers;
  hasDots = isfield(images, 'dotIm');

  if hasDots
    dotIm = images.dotIm;
  else
      bgCutOut2 = [];
  end

  % if isfield(sets,'opticsFile')
  %     opticsDir = dir(sets.opticsFile);
  %     folder = [opticsDir(1).folder,'/'];
  % else
  %     folder = sets.targetFolder;
  % end

  foundOptics = 0;
  optics = struct();

  if isfield(sets, 'opticsFile') && (isfile(sets.opticsFile) || isfile(fullfile(sets.targetFolder, sets.opticsFile)))

    try
      [optics.NA, optics.pixelSize, optics.waveLength] = get_optic_params(fullfile(sets.targetFolder, sets.opticsFile));
    catch
      [optics.NA, optics.pixelSize, optics.waveLength] = get_optic_params(sets.opticsFile);
    end

    optics.sigma = 1.22 * optics.waveLength / (2 * optics.NA) / optics.pixelSize; % Calculate width of PSF
    optics.logSigma = optics.sigma;
    foundOptics = 1;
  else
    optics.NA = nan;
    optics.waveLength = nan;
  end

  if isfield(sets, 'logSigmaNm') && isfield(sets, 'pxnm')

    if not(foundOptics)
      optics.sigma = sets.logSigmaNm / sets.pxnm;
    end

    optics.logSigma = sets.logSigmaNm / sets.pxnm;
    optics.pixelSize = sets.pxnm;
    foundOptics = 1;
  end

  if not(foundOptics)
    throw(MException('segment:optics', 'No optics file or optics settings found'))
  end

  % if imageNumber == 1
  %     fprintf('Width of point spread function estimated to be %.2f pixels.\n',optics.sigma);
  % end

  % Use _some_ scoring to segment the image as a BW image with molecules in white
  %%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%
  % edgePx = 3; %Arbitrary minimum distance to edge (filters away regions that are less than edgePx pixels from the edge of the image)
  %%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%

  tic;
  % Filter image with LoG filter
  n = ceil(6 * optics.logSigma);
  n = n + 1 -mod(n, 2);
  filt = fspecial('log', n, optics.logSigma);

  if ~isa(imAverage, 'double')
    imAverage = double(imAverage);
    fprintf('Averaged image has been converted to type "double"\n');
  end

  logim = imfilter(imAverage, filt);

  % Find zero crossing contours
  if sets.showMolecules
%     figure(1 + (imageNumber - 1) * 5)
%     t = tiledlayout(hPanelResult,4,2,'TileSpacing','tight','Padding','tight');
    axes(tiles.logFilt);
    imshow(logim, 'InitialMagnification', 'fit')
 
%     ax1 = uiaxes(hPanelResult);
%     imshow(logim, 'InitialMagnification', 'fit','Parent',ax1)
    title([imageName, ' LoG filtered'])
  end

    thedges = imbinarize(logim, 0);
    thedges(1:end,[ 1 end]) = 1; % things around the boundary should also be considered
    thedges([ 1 end],1:end) = 1;

  [B, L] = bwboundaries(thedges, 'holes');
  
  % assume background is the biggest detected feature, then order does not
  % matter. TODO:If not true, introduce a check for consistency with respect to
  % intensity. Check 1 max, 2 max, 3rd max
  
  % make this quicker?
%   featureSizes = arrayfun(@(x) sum(L==x,'all'),1:length(B));
%     featureSizes = arrayfun(@(x) size(B,1),1:length(B));
  featureSizes = arrayfun(@(x) sum(L==x,'all'),0:10);

  [a,b] = sort(featureSizes,'desc');
  % mean of top sizes
%   [minMean,pos] = min(arrayfun(@(x) mean(imAverage(L==x)),b(1:min(5,end))));
  
  % first feature is usually background, but in case 
%   B(1)=[];  %remove first
  L(L==b(1)-1) = 0;
%   B(b(1)-1) = [];
%   L(1)=[];

  % Calculate edge score around all edges
  [~, Gdir] = imgradient(logim);
  dist = wd_calc(optics.sigma); % Very simple function - perhaps do elsewhere.
  meh = zeros(1, length(B));
  stat = @(h) mean(h); % This should perhaps be given from the outside
  fprintf('Total number of regions: %i.\n', length(B));

  for k = 1:length(B)% Filter out any regions with artifacts in them

    if not(isempty(centers))
      in = inpolygon(centers(:, 1), centers(:, 2), B{k}(:, 2), B{k}(:, 1));
    else
      in = 0;
    end

    if ~in
      meh(k) = edge_score(B{k}, logim, Gdir, dist, stat); %/optics.sigma^3;% What is the point of this division? //Erik
    else
      meh(k) = -Inf;
    end

  end

  % If showScores, show histogram of log of region scores
  if sets.showScores
%     figure(2 + (imageNumber - 1) * 5)
    axes(tiles.molScores);
    h1 = histogram(log(meh(meh > 1)));
    title([imageName, ' edge scores'])
  end

  % If showMolecules, show image to outline molecules
  if sets.showMolecules
    molFigNum = 3 + (imageNumber - 1) * 5;
%     figure(molFigNum)
%     movies.c = nexttile(t);
    axes(tiles.molDet);
    imshow(mat2gray(imAverage), 'InitialMagnification', 'fit');
    title([imageName, ' detected molecules'])

    if hasDots
      dotFigNum = 4 + (imageNumber - 1) * 5;
%       figure(dotFigNum)
      axes(tiles.dotDet);
      imshow(mat2gray(dotIm), 'InitialMagnification', 'fit');
      title([imageName, ' dot image']);
    end

  end

  %%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%
  % lowLim = exp(-4.5);     % Set the low score threshold to consider a region "signal" (very important)
  % highLim = exp(24);   % Arbitrary higher bound not utilized at this point. (ignore this for now)
  % elim = .8;			 % Set lower limit for eccentricity of region (removes dots and circles and keeps long shapes)
  % ratlim = .5;		 % Set lower limit for ratio of area of region to the convex region formed around (removes "wiggly" regions)
  % lengthLims = [50 1000]; % Set lower and upper limit for the length of the molecule (pixels)
  % widthLims = [0 50]; 		% Set lower and upper limit for the width of the molecule (pixels)
  %%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%

  lowLim = sets.lowLim;
  highLim = sets.highLim;
  elim = sets.elim;
  ratlim = sets.ratlim;
  lengthLims = sets.lengthLims;
  widthLims = sets.widthLims;
  sigmaBgLim = sets.sigmaBgLim;

  if sets.autoThreshBars
    logEdgeScores = log(meh(meh > 1));
    lowestScore = min(logEdgeScores);
    [autoThreshRel, em] = graythresh(logEdgeScores(:) - lowestScore);
    autoThresh = exp(autoThreshRel * (max(logEdgeScores) - lowestScore) + lowestScore);
    lowLim = autoThresh;
  end

  if sets.showScores
%     figure(2 + (imageNumber - 1) * 5)
    axes(tiles.molScores)
    hold on
    %    line([log(autoThresh) log(autoThresh)],[0 max(h1.Values)],'LineStyle','--','Color','red','LineWidth',2)
    if sets.autoThreshBars
      mess = 'Automated threshold';
      line([log(lowLim) log(lowLim)], [0 max(h1.Values)], 'LineStyle', '--', 'Color', 'red', 'LineWidth', 2)
    else
      line([log(lowLim) log(lowLim)], [0 max(h1.Values)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 2)
      mess = 'User set threshold';
    end

    text(1.1 * log(lowLim), 2/3 * max(h1.Values), mess, 'FontSize', 14)
    hold off
  end

  % Filter molecules via the tweak parameters above
  accepted = 0;
  D = B;
  newL = zeros(size(L));
  trueedge = cell(1, length(B));
  scores = nan(1, length(B));

  % Background stuff
  bgPixels = imAverage(L ==0);    
  bgMean = trimmean(bgPixels(:), 10);
  
    if length(bgPixels) > 10000
        bgPixels(bgPixels>bgMean+2*std(bgPixels)) = [];
    end

  bgCutOut = nan(size(imAverage));
  bgCutOut(L ==0) = imAverage(L ==0);
  
%     if sets.showMolecules
%         %     figure(1 + (imageNumber - 1) * 5)
%         %     t = tiledlayout(hPanelResult,4,2,'TileSpacing','tight','Padding','tight');
%         axes(tiles.bg);
%         imagesc(bgCutOut)
%         colormap(gray)
%         %     ax1 = uiaxes(hPanelResult);
%         %     imshow(logim, 'InitialMagnification', 'fit','Parent',ax1)
%         title([imageName, ' bg pixels'])
%     end
% 


  if sigmaBgLim > 0
    bgStd = sqrt(trimmean(bgPixels(:).^2, 10) - bgMean.^2);
  end

  bgSubtractedIm = cellfun(@(x) x- bgMean,registeredIm,'un',false);
  for ix=1:length(bgSubtractedIm)
      bgSubtractedIm{ix}(bgSubtractedIm{ix} <= 0) = nan;
  end
  
    if sets.showMolecules
        axes(tiles.molDet);
        hold on
    end
        % what to do with the first??
for k = 1:length(B)% Filter any edges with lower scores than lim
    acc = mol_filt(B{k}, meh(k), lowLim, highLim, elim, ratlim, lengthLims, widthLims);

    if sigmaBgLim > 0
      numPixelsInMol = sum(L == k, 'all');
      intensAcc = sum(imAverage(L == k)) >= numPixelsInMol * bgMean + sqrt(numPixelsInMol) * sigmaBgLim * bgStd;
    else
      intensAcc = 1;
    end

    if acc && intensAcc
      accepted = accepted + 1;
      trueedge{accepted} = D{k};
      scores(k) = meh(k);
      newL(L == k) = accepted;

      if sets.showMolecules
%         figure(molFigNum)

        plot(trueedge{accepted}(:, 2), trueedge{accepted}(:, 1));
      end

    else
      D{k} = [];
    end

end
  
      if sets.showMolecules
          hold off
      end


  trueedge(accepted + 1:end) = [];
  D = D(~cellfun('isempty', D)); % Remove empty entries (where molecules have been filtered)
  scores = scores(~isnan(scores));

  if hasDots && sets.showMolecules
    axes(tiles.dotDet);
    hold on

    for k = 1:length(trueedge)
%       figure(dotFigNum)
      plot(trueedge{k}(:, 2), trueedge{k}(:, 1));
    end
    
    hold off

  end

  fprintf('Found %i molecules with log edge score >%.2f, eccentricity >%.2f, aRat >%.2f, %i< length <%i, and %i< width <%i.\n', accepted, log(lowLim), elim, ratlim, lengthLims(1), lengthLims(2), widthLims(1), widthLims(2));
  t = toc;
  fprintf('Segmentation completed in %.1f seconds.\n', t);

  % Store each molecules in its own image (this might be very slow)
  folderName = subsref(dir(sets.targetFolder), substruct('.', 'folder'));

  if hasDots
    [molM, bwM, dotM, pos] = generate_molecule_images_fast(D, newL, bgSubtractedIm, dotIm, folderName, runNo, imageName, sets.edgeMargin, sets);
  else
    [molM, bwM, dotM, pos] = generate_molecule_images_fast(D, newL, bgSubtractedIm, [], folderName, runNo, imageName, sets.edgeMargin, sets);
  end

  movies.imageName = imageName;
  movies.molM = molM;
  movies.bwM = bwM;
  movies.pos = pos;

  if sets.showMolecules
    movies.molFigNum = molFigNum;
  end

  if hasDots
    movies.dotFigNum = dotFigNum;
    movies.dotM = dotM;
    bgCutOut2 = nan(size(dotIm)); % background
    bgCutOut2(L ==0) = dotIm(L ==0);

    if sets.showMolecules
    %     figure(1 + (imageNumber - 1) * 5)
    %     t = tiledlayout(hPanelResult,4,2,'TileSpacing','tight','Padding','tight');
    axes(tiles.bg);
    imagesc(bgCutOut2);
    colormap(gray);
    %     ax1 = uiaxes(hPanelResult);
    %     imshow(logim, 'InitialMagnification', 'fit','Parent',ax1)
    title([imageName, ' bg pixels'])
    end

    
    bgMean2 =  trimmean(bgCutOut2(:), 10);
    if length(bgCutOut2(:)) > 10000
        bgCutOut2 = bgCutOut2(bgCutOut2<bgMean2+2*nanstd(bgCutOut2(:)));
    end

  end
% 
%   
%     % Background stuff
%   bgPixels = imAverage(L ==0);    
%   bgMean = trimmean(bgPixels(:), 10);
%   
% 
%   bgCutOut = nan(size(imAverage));
%   bgCutOut(L ==0) = imAverage(L ==0);
%   
% 
% 
