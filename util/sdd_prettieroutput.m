function [ molecules ] = sdd_prettieroutput(output, label )

% This puts the output info in a better way

% Initialize output structure
molecules = struct( ...
    'ImageName', {}, ...
    'ImageID', {}, ...
    'Mol', {}, ...
    'MolID', {}, ...
    'rawBarcode', {}, ...
    'length', {}, ...
    'averageIntensity', {}, ...
    'stdev',{},...
    'continuous', {});

% Global molecule counter
MolID = 0;

% Loop through images
for i = 1:numel(output)

    % Image information
    imageName = output{i}.name;
    if nargin <2
    label=imageName(1:3);
    end
    expBars = output{i}.expBars;

    % Loop through molecules in the image
    for m = 1:numel(expBars)

        MolID = MolID + 1;

        % Extract raw intensity profile
        rawBarcode = expBars{m}.rawBarcode;

        % Basic molecule properties
        N = numel(rawBarcode);
        averageIntensity = mean(rawBarcode, 'omitnan');

        % Check whether the barcode contains internal NaNs
        % (NaNs at the beginning or end are allowed)
        validIdx = find(~isnan(rawBarcode));

        if isempty(validIdx)
            continuous = 0;
        else
            continuous = ~any(isnan(rawBarcode(validIdx(1):validIdx(end))));
        end

        % Store molecule information
        molecules(MolID).ImageName = imageName;
        molecules(MolID).ImageID = i;
        molecules(MolID).Mol = m;
        molecules(MolID).MolID = MolID;
        molecules(MolID).rawBarcode = rawBarcode;
        molecules(MolID).length = N;
        molecules(MolID).averageIntensity = averageIntensity;
        molecules(MolID).averageIntensity = nanstd(rawBarcode);
        molecules(MolID).continuous = continuous;

        molecules(MolID).treatment=label;

    end

end
end


