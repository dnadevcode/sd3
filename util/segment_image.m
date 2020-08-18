function [movies,scores,optics,lengthLims,lowLim,widthLims] = segment_image(images,imageName,imageNumber,runNo,experiment,actions,sets)

registeredIm = images.registeredIm;
imAverage = images.imAverage;
%imDenoised = images.imDenoised;
centers = images.centers;
hasDots = isfield(images,'dotIm');
if hasDots
	dotIm = images.dotIm;
end
% if isfield(experiment,'opticsFile')
%     opticsDir = dir(experiment.opticsFile);
%     folder = [opticsDir(1).folder,'/'];
% else
%     folder = experiment.targetFolder;
% end

foundOptics = 0;
optics = struct();
if isfield(experiment, 'opticsFile') && (isfile(experiment.opticsFile) || isfile(fullfile(experiment.targetFolder, experiment.opticsFile)))
    try
        [optics.NA,optics.pixelSize,optics.waveLength] = get_optic_params(fullfile(experiment.targetFolder, experiment.opticsFile));
    catch
        [optics.NA,optics.pixelSize,optics.waveLength] = get_optic_params(experiment.opticsFile);
    end
    optics.sigma = 1.22*optics.waveLength/(2*optics.NA) / optics.pixelSize; % Calculate width of PSF
    if imageNumber == 1
        fprintf('Width of point spread function estimated to be %.2f pixels.\n',optics.sigma);
    end
    foundOptics = 1;
else
    optics.NA = nan;
    optics.waveLength = nan;
end
if isfield(experiment, 'psfnm') && isfield(experiment, 'pxnm')
    optics.sigma = experiment.psfnm/experiment.pxnm;
    optics.pixelSize = experiment.pxnm;
    foundOptics = 1;
end
if not(foundOptics)    
    throw(MException('segment:optics', 'No optics file or optics settings found'))
end

% Use _some_ scoring to segment the image as a BW image with molecules in white
%%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%
% edgePx = 3; %Arbitrary minimum distance to edge (filters away regions that are less than edgePx pixels from the edge of the image) 
%%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%

tic;
% Filter image with LoG filter
n = ceil(6*optics.sigma);
n = n + 1 -mod(n,2);
filt = fspecial('log',n,optics.sigma);
if ~isa(imAverage,'double')
	imAverage = double(imAverage);
	fprintf('Averaged image has been converted to type "double"\n');
end
logim = imfilter(imAverage,filt);  

% Find zero crossing contours
if actions.showMolecules
	figure(1+(imageNumber-1)*5)
	imshow(logim,'InitialMagnification','fit') 
	title([imageName,' LoG filtered'])
end

thedges = imbinarize(logim,0);
[B,L] = bwboundaries(thedges,'holes');

% Calculate edge score around all edges
[~,Gdir] = imgradient(logim);
dist = wd_calc(optics.sigma); % Very simple function - perhaps do elsewhere.
meh = zeros(1,length(B));
stat = @(h) mean(h);  % This should perhaps be given from the outside
fprintf('Total number of regions: %i.\n',length(B));
for k = 1:length(B) % Filter out any regions with artifacts in them
    if not(isempty(centers))
        in = inpolygon(centers(:,1),centers(:,2),B{k}(:,2),B{k}(:,1));
    else
        in = 0;
    end
    if ~in
        meh(k) = edge_score(B{k},logim,Gdir,dist,stat);%/optics.sigma^3; % What is the point of this division? //Erik
    else
        meh(k) = -Inf;
    end
end

% If showScores, show histogram of log of region scores
if actions.showScores
	figure(2+(imageNumber-1)*5)
	h1 = histogram(log(meh(meh>0)));
	title([imageName,' edge scores'])
end

% If showMolecules, show image to outline molecules
if actions.showMolecules
    molFigNum = 3+(imageNumber-1)*5;
    figure(molFigNum)
    imshow(mat2gray(imAverage),'InitialMagnification','fit');
	title([imageName,' detected molecules'])
	if hasDots
		dotFigNum = 4+(imageNumber-1)*5;
		figure(dotFigNum)
		imshow(mat2gray(dotIm),'InitialMagnification','fit');
		title([imageName,' dot image']);
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

lowLim = experiment.lowLim;
highLim = experiment.highLim;
elim = experiment.elim;
ratlim = experiment.ratlim;
lengthLims = experiment.lengthLims;
widthLims = experiment.widthLims;
sigmaBgLim = experiment.sigmaBgLim;

if actions.autoThreshBars
    logEdgeScores = log(meh(meh>0));
    lowestScore = min(logEdgeScores);
    [autoThreshRel, em] = graythresh(logEdgeScores-lowestScore);
    autoThresh = exp(autoThreshRel*(max(logEdgeScores)-lowestScore)+lowestScore);
	lowLim = autoThresh;
end

if actions.showScores
    figure(2+(imageNumber-1)*5)
    hold on
%    line([log(autoThresh) log(autoThresh)],[0 max(h1.Values)],'LineStyle','--','Color','red','LineWidth',2)
	if actions.autoThreshBars
		mess = 'Automated threshold';
    	line([log(lowLim) log(lowLim)],[0 max(h1.Values)],'LineStyle','--','Color','red','LineWidth',2)
	else
    	line([log(lowLim) log(lowLim)],[0 max(h1.Values)],'LineStyle','--','Color','black','LineWidth',2)
		mess = 'User ser threshold';
	end
    text(1.1*log(lowLim),2/3*max(h1.Values),mess,'FontSize',14)
    hold off
end

% Filter molecules via the tweak parameters above
accepted = 0;
D = B;
newL = zeros(size(L));
trueedge = cell(1,length(B));
scores = nan(1,length(B));

if sigmaBgLim > 0
    bgPixels = imAverage(L == 0);
    bgMean = trimmean(bgPixels(:), 10);
    bgStd = sqrt(trimmean(bgPixels(:).^2, 10)-bgMean.^2);
    bgThresh = bgMean + sigmaBgLim*bgStd;
end

for k = 1:length(B) % Filter any edges with lower scores than lim
	acc = mol_filt(B{k},meh(k),lowLim,highLim,elim,ratlim,lengthLims,widthLims);
    if sigmaBgLim > 0
        numPixelsInMol = sum(L == k, 'all');
        intensAcc = sum(imAverage(L == k)) >= numPixelsInMol*bgThresh;
    else
        intensAcc = 1;
    end
	if acc && intensAcc
		accepted = accepted + 1;
	    trueedge{accepted} = D{k};
		scores(k) = meh(k);
		newL(L == k) = accepted;
        if actions.showMolecules
			figure(molFigNum)
            hold on
            plot(trueedge{accepted}(:,2),trueedge{accepted}(:,1));
            hold off
        end
	else		
		D{k} = [];
	end
end
trueedge(accepted+1:end) = [];
D = D(~cellfun('isempty',D)); % Remove empty entries (where molecules have been filtered)
scores = scores(~isnan(scores));
if hasDots && actions.showMolecules
	for k = 1:length(trueedge)
		figure(dotFigNum)
    	hold on
	    plot(trueedge{k}(:,2),trueedge{k}(:,1));
    	hold off
	end
end
fprintf('Found %i molecules with log edge score >%.2f, eccentricity >%.2f, aRat >%.2f, %i< length <%i, and %i< width <%i.\n',accepted,log(lowLim),elim,ratlim,lengthLims(1),lengthLims(2),widthLims(1),widthLims(2));
t = toc;
fprintf('Segmentation completed in %.1f seconds.\n',t);

% Store each molecules in its own image (this might be very slow)
folderName = subsref(dir(experiment.targetFolder), substruct('.', 'folder'));
if hasDots
	[molM,bwM,dotM,pos] = generate_molecule_images_fast(D,newL,registeredIm,dotIm,folderName,runNo,imageName,sets.edgeMargin,actions);
else
	[molM,bwM,dotM,pos] = generate_molecule_images_fast(D,newL,registeredIm,[],folderName,runNo,imageName,sets.edgeMargin,actions);
end
movies.imageName = imageName;
movies.molM = molM;
movies.bwM = bwM;
movies.pos = pos;
if actions.showMolecules
    movies.molFigNum = molFigNum;
end
if hasDots 
	movies.dotM = dotM;
    if actions.showMolecules
        movies.dotFigNum = dotFigNum;
    end
end
