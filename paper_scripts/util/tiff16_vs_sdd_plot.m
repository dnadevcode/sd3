function [B, out1,hPanelResult1,images1,movies1,barcodes1] = tiff16_vs_sdd_plot(foldRun1,foldRun,b,keepRows,idx)

%%
% idx = 4; % good
B = imread(fullfile(foldRun1, ['RG (', num2str(b(keepRows(idx))), ').jpeg']));


dFSingle = {fullfile(foldRun,['Green (' num2str(b(keepRows(idx))), ').tiff']) };
[out1,hPanelResult1,images1,movies1,barcodes1]  = sdd_script('sdd_fig4.txt',[],dFSingle{1});




end

