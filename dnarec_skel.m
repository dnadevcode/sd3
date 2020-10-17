function output = dnarec_skel(target, inputsets)

% Store data path and assumptions on experimental conditions

experiment = struct();
actions = struct();
functions = struct();
sets = struct();

if nargin > 0
    experiment.targetFolder = target;
else
    experiment.targetFolder = pwd;
end
if nargin < 2
    
    %%% User default settings %%%
    experiment.barFlag = ''; % Flag for specifying that the image is an image of barcodes. Only images containing the flag in their name will be analysed
    experiment.dotFlag = 'C=0'; % Flag for specifying that the image is a dot image. The flag must be contanied in the image name.
    
    experiment.lowLim = exp(0); % Set the low score threshold to consider a region "signal" (very important)
    experiment.elim = .8; % Set lower limit for eccentricity of region (removes dots and circles and keeps long shapes)
    experiment.ratlim = .4; % Set lower limit for ratio of area of region to the convex region formed around (removes "wiggly" regions)
    experiment.lengthLims = [50 inf]; % Set lower and upper limit for the length of the molecule (pixels)
    experiment.widthLims = [1 inf]; % Set lower and upper limit for the width of the molecule (pixels)
    experiment.dotScoreMin = 1e4; % Dot score lower threshold.
    
    actions.showScores = 0; % Show a histogram of the scores for the regions as well as the dots, which helps setting tweak parameters
    actions.showMolecules = 0; % Show plots of detected molecules and images
    actions.saveMolecules = 0; % Save individual molecule and dot images (2D)
    actions.saveBars = 0; % Save the generated barcodes and dot barcodes
    actions.autoThreshBars = 0;
    actions.autoThreshDots = 0;
    %%% End of settings %%%
    
else
    experiment.barFlag = inputsets.barFlag;
    experiment.dotFlag = inputsets.dotFlag;
    
    experiment.logSigmaNm = inputsets.logSigmaNm;
    experiment.pxnm = inputsets.pxnm;
    
    experiment.lowLim = inputsets.lowLim;
    experiment.elim = inputsets.elim;
    experiment.ratlim = inputsets.ratlim;
    experiment.lengthLims = inputsets.lengthLims;
    experiment.widthLims = inputsets.widthLims;
    experiment.dotScoreMin = inputsets.dotScoreMin;
    
    actions.showScores = inputsets.showScores;
    actions.showMolecules = inputsets.showMolecules;
    actions.saveMolecules = inputsets.saveMolecules;
    actions.saveBars = inputsets.saveBars;
    actions.autoThreshBars = inputsets.autoThreshBars;
    actions.autoThreshDots = inputsets.autoThreshDots;
end

% Legacy settings (remove?)
experiment.opticsFile = 'optics.txt'; % Optionally specified optics file location (full path needed). If this is specified, the software assumes the optical setup to be identical for all images
experiment.highLim = inf; % Arbitrary higher bound not utilized at this point. (ignore this for now)
sets.edgeMargin = 3; % Minimum distance (in pixels) from image edge for a molecule to be included in the analysis.
sets.dotMargin = 1; % Minimum distance (in pixels) from molecule end for a dot to be included in the analysis.
sets.deltaCut = 1; % Number of sigma_psf uncertainty for extract_barcodes.
sets.fragLengthRangeBp = [4 8 12]; % Specfiy range breakpoints (micrometers), for the number of DNA fragments in each range.

% Experimental settings
experiment.sigmaBgLim = 0; % Set lower limit for number of standard deviations in molecule intensity from the background


% targetFile = ''; % Optional specified target image. If the whole folder should be analysed, leave "targetFile = '';"

%%% IFF only one (DNA molecule) image should be analyzed, write the full filename here
%targetFile = 'Snap-1283.czi - C=1.tif';
%%% else, comment the line.

% Specify functions to be used for analysis. Set these equal 1 to have the action performed or 0 to skip the action.
% You do not need to worry about these unless you would like to avoid showing or saving images / barcodes.

%  actions.denoiseImage = 0; 			% Denoise the image using trend removal in 2D (shouldnt be needed for benchtop microscopy)
%  actions.removeRegions = 0;			% Remove regions outside a central brigth region (can remove too much if illumination is uniform)
%  actions.checkSameOrientation = 1;	% Check that images from the same stack (of the same subject) are identically oriented. Only applies or image stacks
%  actions.makeSameSize = 1;			% Makes images in imagestacks the same size
% actions.computeFreeConcentrations = 1;
%  actions.getConsensus = 1;			%


%%%% Note for users %%%%
% The rest of the parameters here are only used for automatically running
% p-value calculation and the consensus routine. They can be ignored if only generating barcode .tifs and movies of the corresponding molecules

%  experiment.concY = 4e-2; 		% YOYO1 concentration in units of 1E-6 molar
%  experiment.concN = 6; % Netropsin concentration in units of 1E-6 molar
% % exp.concDNA = 0.2; 		% DNA concentration in units of 1E6 molar
% % exp.seqence_path = PATH_OF_TARGET_SEQUENCE % Specify this if a custom sequences is to be analysed separately


% % Specify settings (untrusted regions, model for theoretical barcodes etc.)
%  import SAD.cb_model;
%  sets.model = cb_model();
%  sets.pval.nulModelPath = 'util/';
%  sets.pval.nulModelName = 'meanFFT.mat';
%  sets.pval.numRandBarcodes = 100;
% import SAD.import_sequence;
% sets.dnaSequence = import_sequence('./');
%  sets.isLinear = 1;
%  if isempty(targetFile)
import SAD.import_images_simple
functions.img_import = @() import_images_simple(experiment);
%  else
% 	 import SAD.import_single_image
% 	 functions.img_import = @() import_single_image(targetFile,experiment,actions);
%  end
import SAD.denoise_images;
functions.img_denoise = @(images) denoise_images(images,actions);
functions.img_segment = @(images,imageNames,imageNumber,saveName) segment_image(images,imageNames,imageNumber,saveName,experiment,actions,sets);
import SAD.extract_barcodes;
functions.bc_extract = @(movies,optics,lengthLims,imageNumber,saveName) extract_barcodes(movies,optics,lengthLims,imageNumber,saveName,sets,experiment,actions);
%  import SAD.analyse_barcodes;
%  functions.bc_analyse = @(barcodes,optics,lengthLims) analyse_barcodes(barcodes,optics,lengthLims,sets,experiment);
%  import SAD.generate_consensus;
%  functions.consensus_gen = @(barcodes) generate_consensus(barcodes);

import SAD.dna_process_folder
output = dna_process_folder(experiment,functions,actions,sets);

folderName = subsref(dir(experiment.targetFolder), substruct('.', 'folder'));
savePath = fullfile(folderName,'dnarecoutput');
save(savePath,'output')
save(savePath,'actions','-append')
save(savePath,'experiment','-append')
