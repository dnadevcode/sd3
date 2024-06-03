function [printName,printName2] = csv_print(output, sets, runNo, i, filtered)
    % prints results as csv
    
    if nargin < 5 || isempty(filtered)
        filtered = 0;
    end
    
    printName ='';
    printName2 = '';

      hasDots = isfield(output{i}, 'dots');
      hasDots2 = isfield(output{i}, 'dots2');


      if ~hasDots
          return;
      end
    
      import SAD.run_scv_print;
      if hasDots
        dotsOut = output{i}.dots;
         [printName] = run_scv_print(output,dotsOut, sets, runNo, i, filtered,1);


      end

      if hasDots2
            dotsOut = output{i}.dots2;
         [printName2] = run_scv_print(output,dotsOut, sets, runNo, i, filtered,2);
      end


end


   