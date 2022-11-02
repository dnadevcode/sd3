function tests = test_funs
    tests = functiontests(localfunctions);
end

%results = runtests('exampleTest.m')
% 
% function setupOnce(testCase)  % do not change function name
% % set a new path, for example
%     % create tabs to display results from the analysis
%     hTabgroup = uitabgroup('Parent',hPanelResult{i});
% 
%     hResScores= uitab(hTabgroup, 'title',strcat('Scores'));
%     t = tiledlayout(hResScores,2,2,'TileSpacing','tight','Padding','tight');
%     tiles.molScores = nexttile(t);
%     tiles.dotScores = nexttile(t);
%     tiles.dotScoresFilt = nexttile(t);
%     tiles.bgScores = nexttile(t);
% 
%     hResplot = uitab(hTabgroup, 'title',strcat('Detected molecules'));
%     t = tiledlayout(hResplot,1,2,'TileSpacing','tight','Padding','tight');
%     tiles.molDet = nexttile(t);
%     tiles.dotDet = nexttile(t);
%     linkaxes([ tiles.molDet tiles.dotDet  ])
% 
%     hResFilt = uitab(hTabgroup, 'title',strcat('Filtered'));
%     t = tiledlayout(hResFilt,1,2,'TileSpacing','tight','Padding','tight');
%     tiles.logFilt = nexttile(t);
%     tiles.bg = nexttile(t);
% end
% % 
% % function teardownOnce(testCase)  % do not change function name
% % % change back to original path, for example
% % end


function testFunctionOne(testCase)
	% test markbars
    
    sdd_mark_bars(movies, barcodes,tiles,sets.extractionMethod);

end