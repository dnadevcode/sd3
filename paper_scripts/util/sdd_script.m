function [output,hPanelResult,images,movies,barcodes]  = sdd_script(sddsets,outname,datafold)

if nargin < 2 || isempty(outname)
    outname = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
end

% test SDD_GUI without any GUI.
sets = Core.Default.read_default_sets(sddsets,0);


if nargin >=3
    sets.folder = datafold;
end


    % create tabbed figure
    hFig = figure('Name', ['SDD-dots GUI v'], ...
        'Units', 'normalized', ...
        'OuterPosition', [0 0 1 1], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none', ...
        'ToolBar', 'figure' ...
    );

    hPanel = uipanel('Parent', hFig);
    h = uitabgroup('Parent',hPanel);
    t1 = uitab(h, 'title', 'SDD');
    tsHCC = uitabgroup('Parent',t1);

    [output,hPanelResult,images,movies,barcodes] = sdd_process_folder(sets.folder , sets,tsHCC);
    savePath = fullfile(fileparts(sets.folder), ['dnarecoutput_',outname]);
    save(savePath, 'output','images','movies','barcodes')
    display(savePath)

end

