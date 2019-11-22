function [pvals,maxCorr,ml_params,shifts,flips] = compute_p_values_ind(expData,randLength,sequence,optics,sets,experiment)

nMol = length(expData);
seqL = length(sequence);
pvals = zeros(1,nMol);
maxCorr = zeros(1,nMol);
ml_params = zeros(2,nMol);
shifts = zeros(1,nMol);
flips = zeros(1,nMol);

% Generate random sequences (we assume that there is no inherent length scale in 
% these barcodes except for the psf) 
randSeq = generate_random_sequences(randLength,optics,sets);
nRand = sets.pval.numRandBarcodes;

% for each molecule :
for i = 1:nMol
	% Use the experimental length /bitmask to select partof the random sequences 
   %  (either use original length or change the untrustedPx parameter)
	rawBar = expData{i}.rawBarcode;
	rawL = length(rawBar);
	randBitmask = expData{i}.rawBitmask;
	maxV = zeros(1,nRand);
	% Calculate the theoretical barcode for this length( use interp?)
	[theory.barcode, theory.bitmask,~,~] = compute_cb_theory_ind(sequence,rawL,optics,sets,experiment);
	% for each random sequence
	for j = 1:nRand
		%Calculate correlation coefficients between theory and randoms
		rS = randSeq{j}(1:rawL);
		xcorrs = compute_pcc(rS,theory.barcode,randBitmask,theory.bitmask);
		maxV(j)= max(xcorrs(:));
	end
	% Use maxVals to estimate evd parameters
	params = compute_evd_params(maxV,rawL);
	ml_params(:,i) = [params(1) ; params(2)];
	% Calculate the correlation coefficient between theory and experiment
	xcorrs = compute_pcc(rawBar,theory.barcode,expData{i}.rawBitmask,theory.bitmask);
	[~,tshift,tflip] = get_best_parameters(xcorrs);
	shifts(i) = tshift(1);
	flips(i) = tflip(1);
	maxCorr(i) = max(xcorrs(:));
	% Calculate the pvalue for the molecule (or just extract from the expdata struct)
	if ~isnan(maxCorr(i))
		pvals(i) = 1 - (0.5 + 0.5*(1-betainc(maxCorr(i)^2,0.5,params(1)/2-1,'upper')))^params(2);
	else
		pvals(i) = nan;
	end
end
