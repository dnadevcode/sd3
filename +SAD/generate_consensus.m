function consensusStruct = generate_consensus(barcodes)

	stretchedBarcodes = cellfun(@(x) x.stretchedBarcode,barcodes.expBars,'UniformOutput',false);
	stretchedBitmasks = cellfun(@(x) x.stretchedBitmask,barcodes.expBars,'UniformOutput',false);
	bg = cellfun(@(x) 0,stretchedBarcodes,'UniformOutput',false);
	normset.consensus.barcodeNormalization = 'background';
	tic;
	import CBT.Consensus.Core.make_consensus_as_struct;
	fprintf('Intitating consensus calculations.\n');
   consensusStruct = gen_consensus(stretchedBarcodes,stretchedBitmasks,bg,normset);
%	[consensusStruct,cache] = make_consensus_as_struct(stretchedBarcodes',stretchedBitmasks');
	t = toc;
	fprintf('Consensus barcoding took %.0f seconds.\n',t);

end
