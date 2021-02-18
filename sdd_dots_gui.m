classdef sdd_dots_gui < matlab.apps.AppBase
  
  % Properties that correspond to app components
  properties (Access = public)
    UIFigure matlab.ui.Figure
    RunButton matlab.ui.control.Button
    PathInput matlab.ui.control.EditField
    BrowseButton matlab.ui.control.Button
    
    PxInputLabel matlab.ui.control.Label
    PxInput matlab.ui.control.NumericEditField
    LogSigmaInputLabel matlab.ui.control.Label
    LogSigmaInput matlab.ui.control.NumericEditField
    
    BarFlagLabel matlab.ui.control.Label
    BarFlagInput matlab.ui.control.EditField
    DotFlagLabel matlab.ui.control.Label
    DotFlagInput matlab.ui.control.EditField
    
    EdgeScoreLabel matlab.ui.control.Label
    EdgeScoreInput matlab.ui.control.NumericEditField
    DotScoreLabel matlab.ui.control.Label
    DotScoreInput matlab.ui.control.NumericEditField
    
    WidthMinLabel matlab.ui.control.Label
    WidthMinInput matlab.ui.control.NumericEditField
    WidthMaxLabel matlab.ui.control.Label
    WidthMaxInput matlab.ui.control.NumericEditField
    LengthMinLabel matlab.ui.control.Label
    LengthMinInput matlab.ui.control.NumericEditField
    LengthMaxLabel matlab.ui.control.Label
    LengthMaxInput matlab.ui.control.NumericEditField
    
    EccentricityLabel matlab.ui.control.Label
    EccentricitySlider matlab.ui.control.Slider
    MolRatLabel matlab.ui.control.Label
    MolRatSlider matlab.ui.control.Slider
    
    ScoresCheckBox matlab.ui.control.CheckBox
    ShowMolCheckBox matlab.ui.control.CheckBox
    SaveMolCheckBox matlab.ui.control.CheckBox
    SaveBarCheckBox matlab.ui.control.CheckBox
    AutoBarCheckBox matlab.ui.control.CheckBox
    AutoDotCheckBox matlab.ui.control.CheckBox
  end
  
  properties (Access = private)
    posx = 120;
    posy = 120;
    width = 640;
    height = 480;
    buttonHeight = 30;
    buttonWidth = 90;
    checkbHeight = 25;
    margin = 30;
    paddingW = 30;
    paddingH = 10;
    widthPathInput = 460; % (width-2*margin-buttonWidth-paddingW)
    halfWidth = 275; % (width-2*margin-paddingW)/2
    thirdWidth = 173.33; % (width-2*margin-2*paddingW)/3
    fourthWidth = 122.5; % (width-2*margin-3*paddingW)/4
    
    % Lazy default values for GUI...
    pathInputDefault = pwd;
    pixelInputDefault = 130;
    logSigmaInputDefault = 300;
    barFlagInputDefault = 'C=1';
    dotFlagInputDefault = 'C=0';
    edgeScoreInputDefault = 0;
    dotScoreInputDefault = 10;
    widthMinInputDefault = 1;
    widthMaxInputDefault = inf;
    lengthMinInputDefault = 50;
    lengthMaxInputDefault = inf;
    eccMinInputDefault = 0.8;
    ratMinInputDefault = 0.4;
  end
  
  methods (Access = private)
    
    % Refreshes UIFigure and places it on top
    function RefreshUI(app)
      drawnow;
      app.UIFigure.Visible = 'off';
      app.UIFigure.Visible = 'on';
    end
    
    function AnalysisEnded(app, ~)
      app.RunButton.Enable = true;
      app.RefreshUI();
      msgbox('Analysis complete!');
    end
    
    % Button pushed function: Run Button
    function RunButtonPushed(app, ~)
      app.RunButton.Enable = false;
      app.RefreshUI();
      
      sets = struct();
      sets.pathInput = app.PathInput.Value;
      
      sets.barFlag = app.BarFlagInput.Value;
      sets.dotFlag = app.DotFlagInput.Value;
      
      sets.logSigmaNm = app.LogSigmaInput.Value;
      sets.pxnm = app.PxInput.Value;
      sets.lowLim = exp(app.EdgeScoreInput.Value);
      sets.dotScoreMin = app.DotScoreInput.Value;
      sets.widthLims = [app.WidthMinInput.Value app.WidthMaxInput.Value];
      sets.lengthLims = [app.LengthMinInput.Value app.LengthMaxInput.Value];
      sets.elim = app.EccentricitySlider.Value;
      sets.ratlim = app.MolRatSlider.Value;
      
      sets.showScores = app.ScoresCheckBox.Value;
      sets.showMolecules = app.ShowMolCheckBox.Value;
      sets.saveMolecules = app.SaveMolCheckBox.Value;
      sets.saveBars = app.SaveBarCheckBox.Value;
      sets.autoThreshBars = app.AutoBarCheckBox.Value;
      sets.autoThreshDots = app.AutoDotCheckBox.Value;
      
      save('gui_settings', 'sets');
      
      c = onCleanup(@app.AnalysisEnded);
      
      try
        dnarec_folder_scan(app.PathInput.Value, sets);
      catch ME
        rethrow(ME)
      end
    end
    
    % Button pushed function: Browse Button
    function BrowseButtonPushed(app, ~)
      answer = uigetdir();
      if not(any(answer == 0))
        app.PathInput.Value = answer;
      end
      app.RefreshUI();
    end
  end
  
  % App initialization and construction
  methods (Access = private)
    
    % Create UIFigure and components
    function createComponents(app)
      
      % Create UIFigure
      app.UIFigure = uifigure;
      app.UIFigure.Position = [app.posx app.posy app.width app.height];
      app.UIFigure.Name = 'SDD-dots GUI';
      app.UIFigure.Resize = false;
      
      % Create Run Button
      app.RunButton = uibutton(app.UIFigure, 'push');
      app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
      app.RunButton.Position = [app.width-app.margin-app.buttonWidth ...
        app.margin ...
        app.buttonWidth ...
        app.buttonHeight];
      app.RunButton.Text = 'Run';
      
      % Create PathInput
      app.PathInput = uieditfield(app.UIFigure, 'text');
      app.PathInput.Position = [app.margin ...
        app.height-app.margin-app.buttonHeight ...
        app.widthPathInput ...
        app.buttonHeight];
      app.PathInput.Value = app.pathInputDefault;
      
      % Create Browse Button
      app.BrowseButton = uibutton(app.UIFigure, 'push');
      app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
      app.BrowseButton.Position = [app.width-app.margin-app.buttonWidth ...
        app.height-app.margin-app.buttonHeight ...
        app.buttonWidth ...
        app.buttonHeight];
      app.BrowseButton.Text = 'Browse';
      
      % Create PxInput
      app.PxInputLabel = uilabel(app.UIFigure);
      app.PxInputLabel.Position = [app.margin ...
        app.height-app.margin-app.paddingH-2*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.PxInputLabel.Text = 'Pixel size (nm)';
      app.PxInput = uieditfield(app.UIFigure, 'numeric');
      app.PxInput.Position = [app.margin ...
        app.height-app.margin-app.paddingH-3*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.PxInput.Value = app.pixelInputDefault;
      
      % Create LogSigmaInput
      app.LogSigmaInputLabel = uilabel(app.UIFigure);
      app.LogSigmaInputLabel.Position = [app.margin+app.paddingW+app.thirdWidth ...
        app.height-app.margin-app.paddingH-2*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.LogSigmaInputLabel.Text = 'Width of LoG filter (nm)';
      app.LogSigmaInput = uieditfield(app.UIFigure, 'numeric');
      app.LogSigmaInput.Position = [app.margin+app.paddingW+app.thirdWidth ...
        app.height-app.margin-app.paddingH-3*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.LogSigmaInput.Value = app.logSigmaInputDefault;
      
      % Create BarFlagInput
      app.BarFlagLabel = uilabel(app.UIFigure);
      app.BarFlagLabel.Position = [app.margin ...
        app.height-app.margin-2*app.paddingH-4*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.BarFlagLabel.Text = 'Molecule image flag';
      app.BarFlagInput = uieditfield(app.UIFigure, 'text');
      app.BarFlagInput.Position = [app.margin ...
        app.height-app.margin-2*app.paddingH-5*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.BarFlagInput.Value = app.barFlagInputDefault;
      
      % Create DotFlagInput
      app.DotFlagLabel = uilabel(app.UIFigure);
      app.DotFlagLabel.Position = [app.margin+app.paddingW+app.thirdWidth ...
        app.height-app.margin-2*app.paddingH-4*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.DotFlagLabel.Text = 'Dots image flag';
      app.DotFlagInput = uieditfield(app.UIFigure, 'text');
      app.DotFlagInput.Position = [app.margin+app.paddingW+app.thirdWidth ...
        app.height-app.margin-2*app.paddingH-5*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.DotFlagInput.Value = app.dotFlagInputDefault;
      
      % Create EdgeScoreInput
      app.EdgeScoreLabel = uilabel(app.UIFigure);
      app.EdgeScoreLabel.Position = [app.margin ...
        app.height-app.margin-3*app.paddingH-6*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.EdgeScoreLabel.Text = 'Minimum log(EdgeScore)';
      app.EdgeScoreInput = uieditfield(app.UIFigure, 'numeric');
      app.EdgeScoreInput.Position = [app.margin ...
        app.height-app.margin-3*app.paddingH-7*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.EdgeScoreInput.Value = app.edgeScoreInputDefault;
      
      % Create DotScoreInput
      app.DotScoreLabel = uilabel(app.UIFigure);
      app.DotScoreLabel.Position = [app.margin+app.paddingW+app.thirdWidth ...
        app.height-app.margin-3*app.paddingH-6*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.DotScoreLabel.Text = 'Minimum DotScore';
      app.DotScoreInput = uieditfield(app.UIFigure, 'numeric');
      app.DotScoreInput.Position = [app.margin+app.paddingW+app.thirdWidth ...
        app.height-app.margin-3*app.paddingH-7*app.buttonHeight ...
        app.thirdWidth ...
        app.buttonHeight];
      app.DotScoreInput.Value = app.dotScoreInputDefault;
      
      % Create WidthMinInput
      app.WidthMinLabel = uilabel(app.UIFigure);
      app.WidthMinLabel.Position = [app.margin ...
        app.height-app.margin-4*app.paddingH-8*app.buttonHeight ...
        app.fourthWidth ...
        app.buttonHeight];
      app.WidthMinLabel.Text = 'Minimum width';
      app.WidthMinInput = uieditfield(app.UIFigure, 'numeric');
      app.WidthMinInput.Position = [app.margin ...
        app.height-app.margin-4*app.paddingH-9*app.buttonHeight ...
        app.fourthWidth ...
        app.buttonHeight];
      app.WidthMinInput.Value = app.widthMinInputDefault;
      
      % Create WidthMaxInput
      app.WidthMaxLabel = uilabel(app.UIFigure);
      app.WidthMaxLabel.Position = [app.margin+app.paddingW+app.fourthWidth ...
        app.height-app.margin-4*app.paddingH-8*app.buttonHeight ...
        app.fourthWidth ...
        app.buttonHeight];
      app.WidthMaxLabel.Text = 'Maximum width';
      app.WidthMaxInput = uieditfield(app.UIFigure, 'numeric');
      app.WidthMaxInput.Position = [app.margin+app.paddingW+app.fourthWidth ...
        app.height-app.margin-4*app.paddingH-9*app.buttonHeight ...
        app.fourthWidth ...
        app.buttonHeight];
      app.WidthMaxInput.Value = app.widthMaxInputDefault;
      
      % Create LengthMinInput
      app.LengthMinLabel = uilabel(app.UIFigure);
      app.LengthMinLabel.Position = [app.margin+2*app.paddingW+2*app.fourthWidth ...
        app.height-app.margin-4*app.paddingH-8*app.buttonHeight ...
        app.fourthWidth ...
        app.buttonHeight];
      app.LengthMinLabel.Text = 'Minimum length';
      app.LengthMinInput = uieditfield(app.UIFigure, 'numeric');
      app.LengthMinInput.Position = [app.margin+2*app.paddingW+2*app.fourthWidth ...
        app.height-app.margin-4*app.paddingH-9*app.buttonHeight ...
        app.fourthWidth ...
        app.buttonHeight];
      app.LengthMinInput.Value = app.lengthMinInputDefault;
      
      % Create LengthMaxInput
      app.LengthMaxLabel = uilabel(app.UIFigure);
      app.LengthMaxLabel.Position = [app.margin+3*app.paddingW+3*app.fourthWidth ...
        app.height-app.margin-4*app.paddingH-8*app.buttonHeight ...
        app.fourthWidth ...
        app.buttonHeight];
      app.LengthMaxLabel.Text = 'Maximum length';
      app.LengthMaxInput = uieditfield(app.UIFigure, 'numeric');
      app.LengthMaxInput.Position = [app.margin+3*app.paddingW+3*app.fourthWidth ...
        app.height-app.margin-4*app.paddingH-9*app.buttonHeight ...
        app.fourthWidth ...
        app.buttonHeight];
      app.LengthMaxInput.Value = app.lengthMaxInputDefault;
      
      % Create EccentricitySlider
      app.EccentricityLabel = uilabel(app.UIFigure);
      app.EccentricityLabel.Position = [app.margin ...
        app.height-app.margin-5*app.paddingH-10*app.buttonHeight ...
        app.halfWidth ...
        app.buttonHeight];
      app.EccentricityLabel.Text = 'Minimum eccentricity';
      app.EccentricitySlider = uislider(app.UIFigure);
      app.EccentricitySlider.Position = [app.margin ...
        app.height-app.margin-5*app.paddingH-10*app.buttonHeight ...
        app.halfWidth ...
        3];
      app.EccentricitySlider.Limits = [0 1];
      app.EccentricitySlider.Value = app.eccMinInputDefault;
      
      % Create MolRatSlider
      app.MolRatLabel = uilabel(app.UIFigure);
      app.MolRatLabel.Position = [app.margin+app.paddingW+app.halfWidth ...
        app.height-app.margin-5*app.paddingH-10*app.buttonHeight ...
        app.halfWidth ...
        app.buttonHeight];
      app.MolRatLabel.Text = 'Minimum molecule-to-convex-hull ratio';
      app.MolRatSlider = uislider(app.UIFigure);
      app.MolRatSlider.Position = [app.margin+app.paddingW+app.halfWidth ...
        app.height-app.margin-5*app.paddingH-10*app.buttonHeight ...
        app.halfWidth ...
        3];
      app.MolRatSlider.Limits = [0 1];
      app.MolRatSlider.Value = app.ratMinInputDefault;
      
      % Create ScoresCheckBox
      app.ScoresCheckBox = uicheckbox(app.UIFigure);
      app.ScoresCheckBox.Position = [app.margin+2*app.paddingW+2*app.thirdWidth ...
        app.height-app.margin-app.paddingH-2*app.buttonHeight ...
        app.thirdWidth ...
        app.checkbHeight];
      app.ScoresCheckBox.Text = 'Show score histograms';
      
      % Create ShowMolCheckBox
      app.ShowMolCheckBox = uicheckbox(app.UIFigure);
      app.ShowMolCheckBox.Position = [app.margin+2*app.paddingW+2*app.thirdWidth ...
        app.height-app.margin-2*app.paddingH-2*app.buttonHeight-app.checkbHeight ...
        app.thirdWidth ...
        app.checkbHeight];
      app.ShowMolCheckBox.Text = 'Show detected molecules';
      
      % Create SaveMolCheckBox
      app.SaveMolCheckBox = uicheckbox(app.UIFigure);
      app.SaveMolCheckBox.Position = [app.margin+2*app.paddingW+2*app.thirdWidth ...
        app.height-app.margin-3*app.paddingH-2*app.buttonHeight-2*app.checkbHeight ...
        app.thirdWidth ...
        app.checkbHeight];
      app.SaveMolCheckBox.Text = 'Save detected molecules';
      
      % Create SaveBarCheckBox
      app.SaveBarCheckBox = uicheckbox(app.UIFigure);
      app.SaveBarCheckBox.Position = [app.margin+2*app.paddingW+2*app.thirdWidth ...
        app.height-app.margin-4*app.paddingH-2*app.buttonHeight-3*app.checkbHeight ...
        app.thirdWidth ...
        app.checkbHeight];
      app.SaveBarCheckBox.Text = 'Save barcodes and dots';
      
      % Create AutoBarCheckBox
      app.AutoBarCheckBox = uicheckbox(app.UIFigure);
      app.AutoBarCheckBox.Position = [app.margin+2*app.paddingW+2*app.thirdWidth ...
        app.height-app.margin-5*app.paddingH-2*app.buttonHeight-4*app.checkbHeight ...
        app.thirdWidth ...
        app.checkbHeight];
      app.AutoBarCheckBox.Text = 'Auto-threshold EdgeScore';
      
      % Create AutoDotCheckBox
      app.AutoDotCheckBox = uicheckbox(app.UIFigure);
      app.AutoDotCheckBox.Position = [app.margin+2*app.paddingW+2*app.thirdWidth ...
        app.height-app.margin-6*app.paddingH-2*app.buttonHeight-5*app.checkbHeight ...
        app.thirdWidth ...
        app.checkbHeight];
      app.AutoDotCheckBox.Text = 'Auto-threshold DotScore';
    end
  end
  
  methods (Access = public)
    
    % Construct app
    function app = sdd_dots_gui
      
      if isfile('gui_settings.mat')
        try
          sets = subsref(load('gui_settings'), substruct('.', 'sets'));
          app.pathInputDefault = sets.pathInput;
          app.pixelInputDefault = sets.pxnm;
          app.logSigmaInputDefault = sets.logSigmaNm;
          app.barFlagInputDefault = sets.barFlag;
          app.dotFlagInputDefault = sets.dotFlag;
          app.edgeScoreInputDefault = log(sets.lowLim);
          app.dotScoreInputDefault = sets.dotScoreMin;
          app.widthMinInputDefault = sets.widthLims(1);
          app.widthMaxInputDefault = sets.widthLims(2);
          app.lengthMinInputDefault = sets.lengthLims(1);
          app.lengthMaxInputDefault = sets.lengthLims(2);
          app.eccMinInputDefault = sets.elim;
          app.ratMinInputDefault = sets.ratlim;
        catch ME
          warning(compose('gui_settings.mat could not be read: %s', ...
            ME.message))
        end
      end
      
      % Create and configure components
      createComponents(app)
      
      % Register the app with App Designer
      registerApp(app, app.UIFigure)
      
      if nargout == 0
        clear app
      end
    end
    
    % Code that executes before app deletion
    function delete(app)
      
      % Delete UIFigure when app is deleted
      delete(app.UIFigure)
    end
  end
end