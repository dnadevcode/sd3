function [kymos,lineParams] = fit_line_method(molM,bwM,sPer)
    %   Args:
    %       molM - molecule movies
    %       bwM - molecule bw
    %       sPer - perpendicular parameter, i.e. how many rows to take
    %       perpendicularly

    %   Returns:
    %       kymos - kymographs
    %       lineParams - line parameters

    kymos = cell(1,length(molM));
    lineParams = cell(1,length(molM));

    for i=1:length(molM)
            [k,b] = get_line_parameters(bwM{i});
            kymos{i} = get_kymo(molM{i}, k , b, sPer);
            lineParams{i} = [k b];
    end

    %
end

