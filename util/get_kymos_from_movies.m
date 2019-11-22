function [kymos,lineParams] = get_kymos_from_movies( molM, bwM, sPer)
    % this function gets kymos from movies
    
    kymos = cell(1, length(molM));
    lineParams = cell(1, length(molM));
    
    %import PD.Core.Extraction.get_line_parameters;
    %import PD.Core.Extraction.get_kymo;

    for i=1:length(molM)
         % get line parameters
        [k,b] = get_line_parameters(bwM{i});
        kymos{i} = get_kymo(molM{i}, k , b, sPer);
        lineParams{i} = [k b];
    %    import PD.Core.Extraction.save_kymos;
    %    save_kymos( kymos{i}, fold, i-1  )
    end
end

