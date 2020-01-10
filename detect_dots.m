function dots = detect_dots(dotBars,optics,dotMargin,lengthLims,imageNumber,imageName,actions)

oldsig = optics.sigma;
optics.sigma = 669/509*optics.sigma;
w = 1;
dots = cell(1,numel(dotBars));
% Use optics to define LoG filter
n = ceil(6*optics.sigma);
n = n + 1 - mod(n,2);
range = (n-1)/2;
idx = -range:range;
filt = -1/(pi*optics.sigma^4)*(1-idx.^2/(2*optics.sigma^2)).*exp(-idx.^2/(2*optics.sigma^2));
peaks = cell(1,numel(dotBars));
allscores = [];
sumI = sum(cellfun(@(x) nansum(x),dotBars));
L = sum(cellfun(@(x) sum(numel(x)),dotBars));
averageI = sumI/L;
% Filter images one by one and calculate a score
for i = 1:numel(dotBars)
	if sum(~isnan(dotBars{i})) > lengthLims(1)
        barLength = numel(dotBars{i});
        [initnan,endnan] = nanfind(dotBars{i});
        dotBars{i}(isnan(dotBars{i})) = 0;
		dotBars{i} = dotBars{i}(1+initnan:end-endnan);
		%dotBars{i} = dotBars{i}(~isnan(dotBars{i}));
        logBar = imfilter(dotBars{i},filt);
		[~,peaklocs] = findpeaks(-logBar);
		peaks{i}.scores = zeros(1,numel(peaklocs));
		peaks{i}.locations = peaklocs;
        %peaks{i}.depth = min(peaklocs-initnan-oldsig,barLength-peaklocs-endnan-oldsig);
        peaks{i}.depth = min(peaklocs-oldsig,barLength-peaklocs-endnan-initnan-oldsig); 
        peaks{i}.leftOffset = initnan;
        peaks{i}.rightOffset = endnan;
		for j = 1:numel(peaklocs)
			idx = max(1,round(peaklocs(j)-w*optics.sigma)):min(numel(dotBars{i}),round(peaklocs(j)+w*optics.sigma));
			peaks{i}.scores(j) = sum(dotBars{i}(idx)-averageI);
		end
	allscores = [allscores peaks{i}.scores];
	else
		peaks{i}.scores = [];
		peaks{i}.locations = [];
        peaks{i}.depth = [];
        peaks{i}.leftOffset = [];
        peaks{i}.rightOffset = [];
	end
end

%%%%%%%%%%%%%% TWEAK PARAMETER %%%%%%%%%%%%%%
pmin = 1e4;
%%%%%%%%%%%%%% TWEAK PARAMETER %%%%%%%%%%%%%%
maxScore = max(allscores);
autoThreshRel = graythresh(allscores/maxScore);
autoThresh = maxScore*autoThreshRel;
if actions.autoThreshDots
	pmin = autoThresh;
end
medint = median(allscores(allscores>pmin));
if actions.showScores
    figure(5+(imageNumber-1)*5)
    h2 = histogram(allscores,20);
    title([imageName, 'dot scores'])
    hold on
%    line([autoThresh autoThresh],[0 max(h2.Values)],'LineStyle','--','Color','red','LineWidth',2)
	if actions.autoThreshDots
		mess = 'Automated threshold';
    	line([pmin pmin],[0 max(h2.Values)],'LineStyle','--','Color','red','LineWidth',2)
	else
		mess = 'User set threshold';
    	line([pmin pmin],[0 max(h2.Values)],'LineStyle','--','Color','black','LineWidth',2)
	end
    text(1.1*pmin,2/3*max(h2.Values),mess,'FontSize',14)
    hold off
end

% Locate peak positions in images with "high" scores
totdots = 0;
marginDots = 0;
for i = 1:numel(dotBars)
	mask = peaks{i}.scores > pmin;
    endMask = peaks{i}.depth > dotMargin;
    accMask = and(mask,endMask);
    marginDots = marginDots + sum(mask.*~endMask);
	dots{i}.locations = peaks{i}.locations(accMask);
    dots{i}.depth = peaks{i}.depth(accMask);
	dots{i}.N = sum(accMask);
	totdots = totdots + dots{i}.N;
	dots{i}.val = peaks{i}.scores(accMask)/medint;
    dots{i}.leftOffset = peaks{i}.leftOffset;
    dots{i}.rightOffset = peaks{i}.rightOffset;
end
fprintf('Found %i dots.\n',totdots);
fprintf('Rejected %i dots due to vicinity to molecule end.\n',marginDots);
end

function [initnan,endnan] = nanfind(bar)
    stillNaN = isnan(bar(1));
    i = 0;
    while stillNaN && i < numel(bar)
        i = i + 1;
        stillNaN = isnan(bar(i));
    end
    initnan = i;
    
    i = numel(bar);
    stillNaN = isnan(bar(i));
    while stillNaN && i > 1
        i = i - 1;
        stillNaN = isnan(bar(i));
    end
    endnan = numel(bar)-i;
end
