function shouldAbort = prompt_figure_excess(numBars, numDots, actions)

shouldAbort = false;
numFiguresScoresPerBars = 1; % Number of figures shown per image and action
numFiguresMoleculesPerBars = 2; % Number of figures shown per image and action
numFiguresScoresPerDots = 1; % Number of figures shown per image and action
numFiguresMoleculesPerDots = 1; % Number of figures shown per image and action
threshNumFigures = 10; % If this number of figures will be opened, prompt user to abort

numFiguresPredicted = numBars*(actions.showScores*numFiguresScoresPerBars + actions.showMolecules*numFiguresMoleculesPerBars) + ...
    numDots*(actions.showScores*numFiguresScoresPerDots + actions.showMolecules*numFiguresMoleculesPerDots);

if numFiguresPredicted >= threshNumFigures
    answer = questdlg(compose(strcat("With the current settings %i figures will be shown. ", ...
        "Proceed anyway?"), numFiguresPredicted), ...
        'Excess figures', ...
        'Ok', ...
        'Abort', ...
        'Ok');
    if not(strcmp(answer, 'Ok'))
        shouldAbort = true;
    end
end

end

