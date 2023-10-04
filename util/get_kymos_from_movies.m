function [kymos, lineParams, xyPars,distance] = get_kymos_from_movies( molM, bwM, sPer, method)
    % get_kymos_from_movies this function gets kymos from movies
    
    %   Args:
    %       molM - molecule movie
    %       bwM - molecule bitmask
    %       sPer - how many pixels to take perpendicularly
    %       method - extraction method
    %
    %   Returns:
    %       kymos - extracted kymograph
    %       lineParams - line parameters for extracted kymograph
    %       xyPars - xy parameters for extracted kymo
    
    if nargin < 4
        method = 1;
    end

    % These will be returned as output
    kymos = cell(1, length(molM));
    lineParams = cell(1, length(molM));
    xyPars = cell(1, length(molM));
    distance =  cell(1, length(molM));

    import Core.LineExtraction.fit_line_method;
    import Core.LineExtraction.fit_spline_method;

    %%
    switch method
        case 1 %"linear"
            [kymos,lineParams,xyPars] = fit_line_method(molM,bwM,sPer);
        case 2 % spline
             [kymos,lineParams,xyPars,distance] = fit_spline_method(molM,bwM,sPer);
%         case 2 % match-pairs
        otherwise
            error('No valid alignment method selected');
    end
    
%%
    %   % plot
    plotExample = 0;
    if plotExample ==1
        idx = 1
        %             detailed_analysis_plot(movies,barcodes,idx)
        figure,
        tiledlayout(1,3);nexttile
        imagesc(molM{idx})
        nexttile
        imagesc(bwM{idx})
        hold on
        plot(xyPars{idx}{2},xyPars{idx}{1},'->','MarkerSize',3, 'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
%         nexttile
%         hold on 
%         plot(xF,yF,'redx')      
        colormap(gray)
    end
end

