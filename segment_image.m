function [movies,scores,optics,lengthLims] = segment_image(images,imageName,imageNumber,runNo,experiment,actions)

registeredIm = images.registeredIm;
imAverage = images.imAverage;
%imDenoised = images.imDenoised;
centers = images.centers;
hasDots = isfield(images,'dotIm');
if hasDots
	dotIm = images.dotIm;
end
if isfield(experiment,'opticsFile')
    opticsDir = dir(experiment.opticsFile);
    folder = [opticsDir(1).folder,'/'];
else
    folder = experiment.targetFolder;
end
[optics.NA,optics.pixelSize,optics.waveLength] = get_optic_params(folder);
optics.sigma = 1.22*optics.waveLength/(2*optics.NA) / optics.pixelSize; % Calculate width of PSF
if imageNumber == 1
	fprintf('Width of point spread function estimated to be %.2f pixels.\n',optics.sigma);
end

% Use _some_ scoring to segment the image as a BW image with molecules in white
%%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%
edgePx = 3; %Arbitrary minimum distance to edge (filters away regions that are less than edgePx pixels from the edge of the image) 
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
    if ~isempty(centers) > 0
        in = inpolygon(centers(:,1),centers(:,2),B{k}(:,2),B{k}(:,1));
    else
        in = 0;
    end
    if ~in
        meh(k) = edge_score(B{k},logim,Gdir,dist,stat)/optics.sigma^3;
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
lowLim = exp(6);     % Set the low score threshold to consider a region "signal" (very important)
highLim = exp(24);   % Arbitrary higher bound not utilized at this point. (ignore this for now)
elim = .8;			 % Set lower limit for eccentricity of region (removes dots and circles and keeps long shapes)
ratlim = .4;		 % Set lower limit for ratio of area of region to the convex region formed around (removes "wiggly" regions)
lengthLims = [15 1000]; % Set lower and upper limit for the length of the molecule (pixels)
widthLims = [0 10]; 		% Set lower and upper limit for the width of the molecule (pixels)
%%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%
[autoThreshRel,em] = graythresh(log(meh(meh>0)));
autoThresh = exp(autoThreshRel*max(log(meh(meh>0))));
if actions.autoThreshBars
	lowlim = autoThresh;
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

for k = 1:length(B) % Filter any edges with lower scores than lim
	acc = mol_filt(B{k},meh(k),lowLim,highLim,elim,ratlim,lengthLims,widthLims);
	if acc
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
fprintf('Found %i molecules with "high" edge score, eccentricity >%.2f, aRat >%.2f, %i < length <%i, and width < %i.\n',accepted,elim,ratlim,lengthLims(1),lengthLims(2),widthLims(2));
t = toc;
fprintf('Segmentation completed in %.1f seconds.\n',t);

% Store each molecules in its own image (this might be very slow)
if hasDots
	[molM,bwM,dotM,pos] = generate_molecule_images_fast(D,newL,registeredIm,dotIm,experiment.targetFolder,runNo,imageName,edgePx,actions);
else
	[molM,bwM,dotM,pos] = generate_molecule_images_fast(D,newL,registeredIm,[],experiment.targetFolder,runNo,imageName,edgePx,actions);
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
