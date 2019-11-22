function output = analyse_barcodes(barcodes,optics,sets,exp)

%fprintf('Sequence "%s..." imported with %i total inputs.\n',sets.lambdaSequence(1:4),length(sets.lambdaSequence));
randLength = 2*max(cellfun(@(x) length(x.rawBarcode),barcodes.expBars));
tic;
[output.pvals,~,~,~,~] = compute_p_values_ind(barcodes.expBars,randLength,sets.dnaSequence,optics,sets,exp);
t = toc;
fprintf('Calculated p values in %.1f seconds.\n',t);
