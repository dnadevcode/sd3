function [] = sdd_gui()

    % Simplified GUI for SDD with possibility of updating results based on
    % parameter choice
    
    % create a simpler version of sdd gui that possibly shows output in the
    % same GUI.
        
    %% Add source files to path and initialize settings
    mFilePath = mfilename('fullpath');
    mfolders = split(mFilePath, {'\', '/'});
    [fd, fe] = fileparts(mFilePath);
    % read settings txt
    setsTable  = readtable(fullfile(fd,'sdd_settings.txt'),'Format','%s%s%s');

    outputRes = []; % for selecting good/bad
    savePath = [];

    processFolders = 1; % whether to process single files or folders
    
    warning(''); % empty warning
    addpath(genpath(fullfile(mfolders{1:end - 1})));
    
    [~, lwid] = lastwarn;
    
    if strcmp(lwid, 'MATLAB:mpath:nameNonexistentOrNotADirectory')
        error('Unexpected error when asserting source folder path.')
    end

    %% Generate UI
       mFilePath = mfilename('fullpath');
            mfolders = split(mFilePath, {'\', '/'});
            versionSDD = importdata(fullfile(mfolders{1:end-1},'VERSION'));
        
    % create tabbed figure
    hFig = figure('Name', ['SDD-dots GUI v',versionSDD{1}], ...
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
    hPanelImport = uitab(tsHCC, 'title', 'Dot import tab');

    % Checklist as a loop
    checkItems = setsTable.Var2(14:20);

   % checkbox for things to plot and threshold
    for i = 1:length(checkItems)
        itemsList{i} = uicontrol('Parent', hPanelImport, 'Style', 'checkbox','Value', str2double(setsTable.Var1{13+i}),'String',{checkItems{i}},'Units', 'normal', 'Position', [0.45 .83-0.05*i 0.3 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    end

    checkItems2 = setsTable.Var2(23);
    for i = 1:length(checkItems2)
        itemsList2{i} = uicontrol('Parent', hPanelImport, 'Style', 'checkbox','Value', str2double(setsTable.Var1{23}),'String',{checkItems2{i}},'Units', 'normal', 'Position', [0.6 .83-0.05*i 0.3 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    end
    
    % parameters with initial values
    textItems =  setsTable.Var2(1:11);    
    values = setsTable.Var1(1:11);

    for i=1:6 % these will be in two columns
        positionsText{i} =   [0.2-0.2*mod(i,2) .88-0.1*ceil(i/2) 0.2 0.03];
        positionsBox{i} =   [0.2-0.2*mod(i,2) .83-0.1*ceil(i/2) 0.15 0.05];
    end
    
    for i=7:11 % these will be in two columns
        positionsText{i} =   [0.2*(i-7) .45 0.15 0.03];
        positionsBox{i} =   [0.2*(i-7) .4 0.15 0.05];
    end

    for i=1:length(textItems)
        textListT{i} = uicontrol('Parent', hPanelImport, 'Style', 'text','String',{textItems{i}},'Units', 'normal', 'Position', positionsText{i},'HorizontalAlignment','Left');%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        textList{i} = uicontrol('Parent', hPanelImport, 'Style', 'edit','String',{strip(values{i})},'Units', 'normal', 'Position', positionsBox{i});%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    end
   
    uicontrol('Parent', hPanelImport, 'Style', 'text','String',{'Minimum molecule eccentricity'},'Units', 'normal', 'Position', [0 0.3 0.2 0.03]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    uicontrol('Parent', hPanelImport, 'Style', 'text','String',{'Minimum molecule-to-convex-hull ratio'},'Units', 'normal', 'Position', [0 0.2 0.2 0.03]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]

    sliderList{1} = uicontrol('Parent', hPanelImport, 'Style', 'slider','Value',0.8,'SliderStep',[0.05 0.1],'Units', 'normal', 'Position', [0 0.25 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    sliderList{2} = uicontrol('Parent', hPanelImport, 'Style', 'slider','Value',0.4,'SliderStep',[0.05 0.1],'Units', 'normal', 'Position', [0 0.15 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    sliderValue{1} = uicontrol('Parent', hPanelImport, 'Style', 'edit','String','0.8','Units', 'normal', 'Position', [0.25 0.25 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    sliderValue{2} = uicontrol('Parent', hPanelImport, 'Style', 'edit','String','0.4','Units', 'normal', 'Position', [0.25 0.15 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]

    % sliders
%     https://se.mathworks.com/matlabcentral/answers/499395-display-value-of-a-slider-in-a-text-box-without-guide
    fun = @(~,e)set(sliderValue{1} ,'String',num2str(get(e.AffectedObject,'Value')));
    addlistener(sliderList{1}, 'Value', 'PostSet',fun);
    fun2 = @(~,e)set(sliderValue{2} ,'String',num2str(get(e.AffectedObject,'Value')));
    addlistener(sliderList{2}, 'Value', 'PostSet',fun2);
        
%     funV = @(~,e)set(sliderList{1} ,'Value',str2double(get(e.AffectedObject,'String')));
%     addlistener(sliderValue{1}, 'Value', 'PostSet',funV);
%     funV2 = @(~,e)set(sliderList{2} ,'Value',str2double(get(e.AffectedObject,'String')));
%     addlistener(sliderValue{2}, 'Value', 'PostSet',funV2);

    dotImport = uicontrol('Parent', hPanelImport, 'Style', 'edit','String',fullfile(fileparts(mfilename('fullpath')),'testfolder'),'Units', 'normal', 'Position', [0 0.9 0.5 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    set(dotImport, 'Min', 0, 'Max', 25)% limit to 10 files via gui;

    dotButton = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'Browse folder'},'Callback',@selection,'Units', 'normal', 'Position', [0.6 0.9 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    dotButtonFile = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'Browse file'},'Callback',@selection2,'Units', 'normal', 'Position', [0.7 0.9 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    dotButtonUI = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'uigetfiles'},'Callback',@selection3,'Units', 'normal', 'Position', [0.8 0.9 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    runButton = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'Run'},'Callback',@run,'Units', 'normal', 'Position', [0.7 0.2 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    
    clearButton = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'Clear visual results'},'Callback',@clear_results,'Units', 'normal', 'Position', [0.7 0.1 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]

    selectGood = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'select Good'},'Callback',@select_good,'Units', 'normal', 'Position', [0.7 0.3 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    set(selectGood, 'Enable', 'off');
    %% Browse folder
    
    function selection(src, event)
        [rawNames] = uigetdir(pwd,strcat(['Select folder with file(s) to process']));
        dotImport.String = strcat(rawNames,filesep);
        processFolders = 1;

    end    

    function selection2(src, event)
        [FILENAME, PATHNAME] = uigetfile(fullfile(pwd,'*.*'),strcat(['Select file(s) to process']),'MultiSelect','on');
        dotImport.String = fullfile(PATHNAME,FILENAME);
        if ~iscell(dotImport.String)
            dotImport.String  = {dotImport.String};
        end
        processFolders = 0;

    end    

    function selection3(src, event)
        [FILENAME] = uipickfiles;
        dotImport.String = FILENAME';
        if ~iscell(dotImport.String)
            dotImport.String  = {dotImport.String};
        end
        processFolders = 0;

    end    



    function clear_results(src, event)
        tsHCC.SelectedTab
        child_handles = allchild(tsHCC);
        delete(child_handles(2:end));
    end

    function run(src, event)
        display(['Started analysis sdd_dots v',versionSDD{1}])
        sets.folder = dotImport.String;
        sets.pxnm = str2double(textList{1}.String);
        sets.logSigmaNm =  str2double(textList{2}.String);
        sets.barFlag =   textList{3}.String{1};
        sets.dotFlag =   textList{4}.String{1};

        sets.lowLim = exp(str2double(textList{5}.String));
        
        sets.dotScoreMin = str2double(textList{6}.String);
        sets.widthLims =  [ str2double(textList{7}.String)  str2double(textList{8}.String)];

        sets.lengthLims = [ str2double(textList{9}.String)  str2double(textList{10}.String)];
        sets.dotMargin = str2double(textList{11}.String);

        sets.elim =  str2double(sliderValue{1}.String);
        sets.ratlim =  str2double(sliderValue{2}.String);

        sets.showScores = itemsList{1}.Value;
        sets.showMolecules = itemsList{2}.Value;
        sets.saveMolecules = itemsList{3}.Value;
        sets.saveBars = itemsList{4}.Value;
        sets.autoThreshBars = itemsList{5}.Value;
        sets.autoThreshDots = itemsList{6}.Value; 
        sets.extractionMethod =   itemsList{7}.Value+1; % detects dots on spline
        sets.denoiseDotImages =  itemsList2{1}.Value;
    
        sets.numSigmasAutoThresh =  str2double(setsTable.Var1{21});
        sets.autoThreshDotsMethod =  strtrim(setsTable.Var1{22}); % autothresh method
        sets.lenRandBar =  str2double(setsTable.Var1{23}); % autothresh method

        
        % parameters not set by GUI
        sets.highLim = inf; % Arbitrary higher bound not utilized at this point. (ignore this for now)
        sets.sigmaBgLim = 0; % Set lower limit for number of standard deviations in molecule intensity from the background
        sets.edgeMargin = 3; % Minimum distance (in pixels) from image edge for a molecule to be included in the analysis.
        sets.deltaCut = 1; % Number of sigma_psf uncertainty for extract_barcodes.
        sets.showDotPeaks = 0;
        sets.fragLengthRangeBp = [4 8 12]; % Specfiy range breakpoints (micrometers), for the number of DNA fragments in each range.
        sets.rLims = [12 23]; % lims for circles in an image.
        sets.denoiseImages = 0;
%         sets.denoiseDotImages = 1;
%         sets.denoiseDotImages =  str2double(setsTable.Var1{23}); % whether to denoise dot images

    
        if  processFolders
            dataFolders = search_folder(sets.folder);
        else
            dataFolders = sets.folder;
        end
    
        outputRes = cell(1,numel(dataFolders));
        for i = 1:numel(dataFolders)       
            % MAIN FUNCTION
            [output,hPanelResult] = sdd_process_folder(dataFolders{i}, sets, tsHCC);
            if ~isempty(output)
                outputRes{i} = output;
                folderName = subsref(dir(dataFolders{i}), substruct('.', 'folder'));
                savePath = fullfile(folderName, ['dnarecoutput_',num2str(i)]);
                save(savePath, 'output')
            end

  
  
        end  
        set(selectGood, 'Enable', 'on');

                     
    end

    function select_good(src, event)
        % select good molecules from output / save in filtered results
        % output
        import UI.goodbadtool;
        [outputNew] = goodbadtool([4 1],outputRes{1}.molRunFold,savePath,outputRes{1}.molRunFold);
        
        sets = outputRes{1}.settings;
        
        import SAD.dnarec_print
        resultsName = dnarec_print(outputNew, sets,outputRes{1}.runNo,1);
        
        fprintf('\n-------------------------------------------------------------------\n');
        fprintf('Filtered results saved in %s', resultsName);
        fprintf('\n-------------------------------------------------------------------\n');

        
    end
    
end


%% search folder functions, possibly move to separate function to reduce clutter


% search folder
function dataFolders = search_folder(target)
  possibles{1} = target;
  nDataFolders = 0;

  while numel(possibles) > 0

    for i = 1:numel(possibles)
      isOutput = filter_output(possibles{i});
      [hasSubfolder, subfolders] = filter_master(possibles{i});
      isLink = filter_link(possibles{i});
      isEmpty = not(hasSubfolder) && not(filter_images(possibles{i}));
      %Check if possibles is an output folder
      if isOutput || isLink || isEmpty
        % Remove element from possibles
      else

        if hasSubfolder
          % Add each subfolder to possibles
          possibles = [possibles subfolders];
        end

        nDataFolders = nDataFolders + 1;
        % Add element to datafolders
        dataFolders{nDataFolders} = possibles{i};
      end

      possibles{i} = [];
    end

    % Remove empty entries from possibles
    emptyEntries = cellfun(@(x) isempty(x), possibles);
    possibles(emptyEntries) = [];
  end

end

function isOutput = filter_output(folder)
  isOutput = contains(folder, 'molsandbars') | contains(folder, 'movies');
end

function hasImages = filter_images(folder)
  content = dir(folder);
  fContent = content(not([content.isdir]));
  hasImages = false;

  for f = fContent'

    try
      imfinfo(fullfile(f.folder, f.name));
      hasImages = true;
      return
    catch
      continue
    end

  end

end

function [hasSubfolder, subfolders] = filter_master(folder)
  content = dir(folder);
  subfolders = cell(1, numel(content));
  nSubfolders = 0;
  isFolder = [content.isdir];
  content = content(isFolder);

  for i = 1:numel(content)
    isOutput = filter_output(content(i).name);
    isLink = filter_link(content(i).name);

    if ~(isLink || isOutput)
      nSubfolders = nSubfolders + 1;
      subfolders{nSubfolders} = fullfile(content(i).folder, content(i).name);
    end

  end

  hasSubfolder = nSubfolders > 0;
end

function isLink = filter_link(folder)
  isLink = strcmp(folder, '.') | strcmp(folder, '..');
end

