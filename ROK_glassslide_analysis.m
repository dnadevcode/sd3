%% ROK glass slide analysis
% Use SDD with the following parameters:
% e=0.9
% ratio = 0.27
% length= 60px
% LoG = over 450 (600?)

%% SDD part
% for the following just run SDD for each sample type once
% then load the output for one type and run two lines for each label

label='control';
molecules_control = sdd_prettieroutput(output, label );


label='low';
molecules_low = sdd_prettieroutput(output, label );


label='mid';
molecules_mid = sdd_prettieroutput(output, label );


% Merge all molecules into one structure array
molecules_all = [molecules_control(:); ...
                 molecules_low(:); ...
                 molecules_mid(:)];



% Treatment labels
treatments = {'control', 'low', 'mid'};


% This cuts to all continuous molecules only
molecules_all_og=molecules_all; % we save a copy of the original, maybe useful if you want to count stuff at the end
molecules_all = molecules_all([molecules_all.continuous] == 1);


%% Analysis parameters
% Tweak these to change some of the results

% Peak detection
nSD = 3;                 % Threshold = mean + nSD * SD
minPeakDistance = 5;     % Minimum distance between peaks (pixels) to be count as separate

% Remove edge pixels before analysis
nEdge = 0;               % Pixels removed from each end (if you think is needed)

% Molecule profile visualization
nMoleculesPlot =4;      % How many random molecules per treatment to look at at the end

% Random seed for reproducibility
rng(42); 

%% RUN!


% Analyze intensity peaks in individual molecules


for k = 1:numel(molecules_all)

    % Get raw intensity profile
    barcode = molecules_all(k).rawBarcode;

    % Remove edge pixels
    if numel(barcode) > 2*nEdge
        barcode(1:nEdge) = NaN;
        barcode(end-nEdge+1:end) = NaN;
    else
        barcode(:) = NaN;
    end

    % Identify valid pixels
    valid = ~isnan(barcode) & ~isinf(barcode);

    % Calculate molecule-specific intensity statistics
    profile = barcode(valid);

    if isempty(profile)
        molecules_all(k).PeakCount = 0;
        molecules_all(k).PeakThreshold = NaN;
        molecules_all(k).PeakPositions = [];
        molecules_all(k).PeakIntensities = [];
        continue
    end

    mu = mean(profile);
    sigma = std(profile);

    % Define peak threshold
    threshold = mu + nSD*sigma;

    % Find contiguous regions without NaNs
    edges = diff([false, valid, false]);
    starts = find(edges == 1);
    stops  = find(edges == -1) - 1;

    % Initialize peak properties
    peakLocations = [];
    peakIntensities = [];

    % Analyze each continuous segment independently
    for s = 1:numel(starts)

        segment = barcode(starts(s):stops(s));

        % Find local maxima exceeding threshold
        [pks, locs] = findpeaks(segment, ...
            'MinPeakHeight', threshold, ...
            'MinPeakDistance', minPeakDistance);

        % Convert segment coordinates to molecule coordinates
        locs = locs + starts(s) - 1;

        peakLocations = [peakLocations, locs];
        peakIntensities = [peakIntensities, pks];

    end

    % Normalize peak locations to 0-50% of molecule length
    if numel(peakLocations) > 0

        firstPixel = find(valid, 1, 'first');
        lastPixel  = find(valid, 1, 'last');

        if lastPixel > firstPixel

            relativePosition = ...
                (peakLocations - firstPixel) / (lastPixel - firstPixel);

            peakPositions = 100 * min(relativePosition, ...
                                      1 - relativePosition);

        else
            peakPositions = zeros(size(peakLocations));
        end

    else
        peakPositions = [];
    end

    % Save results
    molecules_all(k).PeakCount = numel(peakLocations);
    molecules_all(k).PeakThreshold = threshold;
    molecules_all(k).PeakPositions = peakPositions;
    molecules_all(k).PeakIntensities = peakIntensities;

end



% Compare molecule characteristics across treatments

% Numeric treatment identifiers
treatmentID = zeros(numel(molecules_all),1);

for t = 1:numel(treatments)
    treatmentID(strcmp(string({molecules_all.treatment}), ...
        treatments{t})) = t;
end

% Extract properties
molLength = [molecules_all.length]';
avgIntensity = [molecules_all.averageIntensity]';

% Treatment colors
treatmentColors = lines(numel(treatments));

% Reproducible jitter
figure('Name','Molecule characteristics','Color','w');

% 1. Molecule length
subplot(1,2,1)
hold on

boxchart(treatmentID, molLength, ...
    'BoxFaceColor',[0.85 0.85 0.85], ...
    'MarkerStyle','none');

for t = 1:numel(treatments)
    idx = treatmentID == t & isfinite(molLength);

    x = t + (rand(sum(idx),1)-0.5)*0.25;

    scatter(x, molLength(idx), 12, treatmentColors(t,:), ...
        'filled', 'MarkerFaceAlpha',0.35, ...
        'MarkerEdgeAlpha',0.25);
end

xticks(1:numel(treatments))
xticklabels(treatments)
ylabel('Molecule length (pixels)')
title('Molecule length')
xlim([0.5 3.5])
box on


% 2. Average intensity
subplot(1,2,2)
hold on

boxchart(treatmentID, avgIntensity, ...
    'BoxFaceColor',[0.85 0.85 0.85], ...
    'MarkerStyle','none');

for t = 1:numel(treatments)
    idx = treatmentID == t & isfinite(avgIntensity);

    x = t + (rand(sum(idx),1)-0.5)*0.25;

    scatter(x, avgIntensity(idx), 12, treatmentColors(t,:), ...
        'filled', 'MarkerFaceAlpha',0.35, ...
        'MarkerEdgeAlpha',0.25);
end

xticks(1:numel(treatments))
xticklabels(treatments)
ylabel('Average intensity')
title('Average intensity')
xlim([0.5 3.5])
box on







% Distribution of peak counts across treatments

peakCount = [molecules_all.PeakCount]';

% Include all peak-count categories, starting at zero
maxPeaks = max(2, max(peakCount));

peakCategories = 0:maxPeaks;

% Count molecules in each category
peakDistribution = zeros(numel(treatments),numel(peakCategories));

for t = 1:numel(treatments)

    idx = treatmentID == t;

    for p = 1:numel(peakCategories)

        peakDistribution(t,p) = ...
            sum(peakCount(idx) == peakCategories(p));

    end

end

% Convert to percentages
peakPercentage = 100 * peakDistribution ./ ...
    sum(peakDistribution,2);

% Plot
figure('Name','Peak abundance','Color','w');

%b = bar(1:numel(peakCategories), peakPercentage, 'stacked');
b = bar(1:numel(treatments), peakPercentage, 'stacked');
hold on

% Color each peak-count category
peakColors = lines(numel(peakCategories));

for p = 1:numel(peakCategories)
    b(p).FaceColor = peakColors(p,:);
end

% Labels
xticks(1:numel(treatments))
xticklabels(treatments)

ylabel('Molecules (%)')
xlabel('Treatment')
title('Distribution of prominent peaks per molecule')

legend(compose('%d peaks',peakCategories), ...
    'Location','eastoutside');

ylim([0 100])
box on






% Compare peak positions across treatments

% Position bins (5% intervals)
positionEdges = 0:5:50;
positionCenters = positionEdges(1:end-1) + 2.5;

% Initialize matrix
% Rows = positional bins, columns = treatments
positionCounts = zeros(numel(positionCenters),3);

for t = 1:3

    idx = treatmentID == t;

    selectedMolecules = molecules_all(idx);

    % Collect all peak positions
    positions = [selectedMolecules.PeakPositions];

    % Count peaks in each positional bin
    counts = histcounts(positions, positionEdges);

    % Normalize by number of molecules, not number of peaks
    positionCounts(:,t) = counts(:) ./ numel(selectedMolecules) * 100;

end

% Plot
figure('Name','Peak position distribution','Color','w');

b = bar(positionCenters, positionCounts, 'grouped');

hold on

for t = 1:3
    b(t).FaceColor = treatmentColors(t,:);
end

xlabel('Position along molecule (%)')
ylabel('Peaks per 100 molecules')
title('Spatial distribution of prominent peaks')

legend(treatments,'Location','northwest')

xticks(0:5:50)
xlim([0 50])

box on%% Plot vertically stacked intensity profiles

% Select molecules containing at least one prominent peak
selectedIdx = find(peakCount > 0);

% Sort by treatment, then number of peaks (descending)
[~,order] = sortrows( ...
    [treatmentID(selectedIdx), -peakCount(selectedIdx)], [1 2]);

selectedIdx = selectedIdx(order);

nProfiles = numel(selectedIdx);






% Compare peak positions across treatments

positionEdges = 0:5:50;
positionCenters = positionEdges(1:end-1) + 2.5;

positionCounts = zeros(numel(positionCenters),3);

for t = 1:numel(treatments)

    idx = treatmentID == t;
    selectedMolecules = molecules_all(idx);

    positions = [selectedMolecules.PeakPositions];

    counts = histcounts(positions, positionEdges);

    positionCounts(:,t) = ...
        counts(:) ./ numel(selectedMolecules) * 100;

end

% Plot
figure('Name','Peak position distribution','Color','w');
hold on

for t = 1:numel(treatments)

    plot(positionCenters, positionCounts(:,t), ...
        '-o', ...
        'Color',treatmentColors(t,:), ...
        'LineWidth',2, ...
        'MarkerSize',6, ...
        'MarkerFaceColor',treatmentColors(t,:));

end

xlabel('Position along molecule (%)')
ylabel('Peaks per 100 molecules')
title('Spatial distribution of prominent peaks')

legend(treatments,'Location','best')

xticks(0:5:50)
xlim([0 50])

grid on
box on






% Select random molecules for profile visualization



selectedIdx = [];

for t = 1:numel(treatments)

    % Molecules from this treatment containing at least one peak
    candidates = find(treatmentID == t & peakCount > 0);

    % Randomly select up to nMoleculesPlot
    nSelect = min(nMoleculesPlot, numel(candidates));

    if nSelect > 0

        randomOrder = randperm(numel(candidates),nSelect);

        selectedIdx = [selectedIdx; candidates(randomOrder)];

    end

end

% Preserve treatment grouping
selectedIdx = selectedIdx(:);
% Plot selected intensity profiles

xGrid = linspace(0,100,401);

verticalSpacing = 5.5;

figure('Name','Selected molecule profiles', ...
    'Color','w', ...
    'Position',[100 100 1200 800]);

hold on

row = 0;

for k = 1:numel(selectedIdx)

    molID = selectedIdx(k);

    % Original intensity profile
    barcode = molecules_all(molID).rawBarcode(:)';

    % Exclude edge pixels
    barcode(1:min(nEdge,numel(barcode))) = NaN;

    if numel(barcode) > nEdge
        barcode(max(1,end-nEdge+1):end) = NaN;
    end

    % Valid pixels
    valid = isfinite(barcode);

    if sum(valid) < 2
        continue
    end

    % Z-score
    mu = mean(barcode(valid));
    sigma = std(barcode(valid));

    if sigma == 0 || ~isfinite(sigma)
        continue
    end

    z = (barcode-mu)/sigma;

    % Limit extreme values for visualization
    z = max(min(z,4),-4);

    % Molecule coordinates
    firstPixel = find(valid,1,'first');
    lastPixel = find(valid,1,'last');

    % Initialize resampled profile
    yGrid = nan(size(xGrid));

    % Identify continuous segments
    edges = diff([false,valid,false]);

    starts = find(edges == 1);
    stops = find(edges == -1)-1;

    for s = 1:numel(starts)

        segmentIdx = starts(s):stops(s);

        xSegment = 100 * ...
            (segmentIdx-firstPixel)/(lastPixel-firstPixel);

        zSegment = z(segmentIdx);

        if numel(segmentIdx) >= 2

            q = xGrid >= xSegment(1) & ...
                xGrid <= xSegment(end);

            yGrid(q) = interp1(xSegment,zSegment, ...
                xGrid(q),'linear');

        else

            [~,q] = min(abs(xGrid-xSegment));

            yGrid(q) = zSegment;

        end

    end

    % Vertical position
    row = row+1;
    yOffset = row*verticalSpacing;

    % Treatment color
    t = treatmentID(molID);

    % Plot profile
    plot(xGrid,yGrid+yOffset, ...
        'Color',treatmentColors(t,:), ...
        'LineWidth',1.2);

    % Annotate original length
    text(103,yOffset, ...
        sprintf('%d px',molecules_all(molID).length), ...
        'FontSize',9, ...
        'VerticalAlignment','middle');

end
% Treatment legend
legendHandles = gobjects(numel(treatments),1);

for t = 1:numel(treatments)
    legendHandles(t) = plot(nan,nan, ...
        'Color',treatmentColors(t,:), ...
        'LineWidth',2);
end

legend(legendHandles,treatments, ...
    'Location','northeast', ...
    'FontSize',10, ...
    'Box','off');
% Formatting
xlabel('Relative molecule position (%)')
ylabel('Individual molecules (vertically offset)')
title('Representative intensity profiles containing prominent peaks')

xticks(0:10:100)
xlim([-2 130])

ylim([0 (row+1)*verticalSpacing])

yticks([])
box on



