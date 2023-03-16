function [] = detailed_analysis_plot(movies,barcodes,idx)

% 
% indexes = {62};
%     for ii=1:length(indexes)
%           idx=find(barcodes.idx==indexes{ii})
%           figure,imagesc(movies.dotM{barcodes.idx(idx)})
%           hold on
%           plot(barcodes.xy{barcodes.idx(idx)}{2},barcodes.xy{barcodes.idx(idx)}{1},'redx')
%           pos = barcodes.dots{idx}.locations+barcodes.nanid(idx);
%           plot(barcodes.xy{idx}{2}(pos),barcodes.xy{idx}{1}(pos),'blackx')
%     %             plot(-barcodes.lineParams{idx}(1)*(1:size(movies.molM{idx},2))+barcodes.lineParams{idx}(2),'redx')
% 
% %           plot(-barcodes.lineParams{idx}(1)*(1:size(movies.molM{idx},2))+barcodes.lineParams{idx}(2),'redx')
%     % %     %   
%         colormap(gray)
% 
% 
% end


    figure,
%     idx = 62;
    mol = movies.molM{idx};
    mol2 = movies.dotM{idx};
    
    xF = barcodes.xy{idx}{1};
    yF = barcodes.xy{idx}{2};
    
    tiledlayout(2,2);nexttile
    imagesc(mol')
    nexttile
    imagesc(mol2')
    hold on
    %             
    plot(xF,yF,'redx')
    % %             
    colormap(gray)
    pos = barcodes.dots{idx}.locations+barcodes.dots{idx}.leftOffset;
%     pos
    plot(barcodes.xy{idx}{1}(pos),barcodes.xy{idx}{2}(pos),'greeno','MarkerSize',20)
      nexttile
    plot(barcodes.expBars{idx}.rawBarcode)
    title('Yoyo barcode')
    xlabel('Position')
    ylabel('Intensity')

    nexttile
    plot(barcodes.dotBars{idx})
    title('Dot barcode')
    xlabel('Position')
    ylabel('Intensity')
