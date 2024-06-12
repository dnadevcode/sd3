% test SDD_GUI without any GUI.
sets = Core.Default.read_default_sets('sdd_fig3.txt',0);

% sets.folder = '/export/scratch/albertas/data_temp/DOTS/031523 h202 data/testFig3/image for slide 7.czi';
% 

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
%     hPanelImport = uitab(tsHCC, 'title', 'Dot import tab');

[output,hPanelResult,images,movies,barcodes] = sdd_process_folder(sets.folder , sets,tsHCC);


f = fig3(images,movies,barcodes,sets)
set(gcf, 'Color', 'w')

% 
%     ft = 'Times';
%     fsz = 10;        
% %     %%%%%%%%%%%%%%%
%     % Your Figure
%     %%%%%%%%%%%%%%%
%     set(findall(gcf,'type','text'), 'FontSize', fsz,'FontName', ft)
%     set(gca,'FontSize', fsz, 'FontName', ft)

print(f, 'Fig3_no_labels.png', '-dpng', '-r400', '-painters');
%% fig 4

% test SDD_GUI without any GUI.
sets = Core.Default.read_default_sets('sdd_fig4.txt',0);

% sets.folder = '/export/scratch/albertas/data_temp/DOTS/031523 h202 data/testFig3/image for slide 7.czi';
% 

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
%     hPanelImport = uitab(tsHCC, 'title', 'Dot import tab');

[output,hPanelResult,images,movies,barcodes] = sdd_process_folder(sets.folder , sets,tsHCC);

% if ~isempty(output)
%     outputRes{jj} = output;
%     folderName = subsref(dir(dataFolders{jj}), substruct('.', 'folder'));
    savePath = fullfile(sets.folder, ['dnarecoutput_',num2str(1)]);
    save(savePath, 'output')
% end
% fig3(images,movies,barcodes)
