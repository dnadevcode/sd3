% Generate synthetic images (via hpfl code) to test the molecule detection and
% dot extraction

rng('default')

[~,~] = mkdir('output')
%% SETTINGS
% CHIP PARS
% for simulation
chipPars.pixelsize = 110;
chipPars.gain = 2;
chipPars.adfactor = 2;
chipPars.offset = 400; % based on nanoimager
chipPars.roNoise = 4;
chipPars.QE = 0.9;
chipPars.c = 0.002;

% SIM PARS
sim.gapBetweenPeaks = 15000;
sim.numberPeaks = 10;
sim.distTwinPeaks = 2000;
sim.stgap = 15000;
sim.nY = 512;
sim.nX = 512;

labelPars.nmbp = 0.23;

noisePars.lambdaBg = 6;  %10          % Poisson parameter for background
noisePars.lambdaSg = 200; %200
noisePars.sigma = 1.5;
noisePars.lambdaYoyo = 1;

%     camera = 'scmos';
camera = 'scmos';

contourLen = 300;
import Simulate.simulate_stretched_contour;
[fluorophorePosPx,fluorophorePosPxYOYO] = simulate_stretched_contour(contourLen,sim.nY/2,sim.numberPeaks);

% grid = 200:5:400;
% fluorophorePosPx = [[grid]' 200*ones(1, length(grid))'];
% fluorophorePosPxYOYO = [[grid]' 200*ones(1, length(grid))'];

timeFrames = 1;
import Simulate.simulate_fluorophores_camera;

finalImageSCMOS = zeros(sim.nY,sim.nX,timeFrames);
finalImageBGYOYO = zeros(sim.nY,sim.nX,timeFrames);
final1D = zeros(timeFrames,sim.nX);
finalImageBGYOYO1D = zeros(timeFrames,sim.nX);
for tf=1:timeFrames
    [finalImageSCMOS(:,:,tf),final1D(tf,:),~,finalImageBG,final1DBG,finalImageBGYOYO(:,:,tf),finalImageBGYOYO1D(tf,:)] = simulate_fluorophores_camera(camera,sim,chipPars,labelPars,noisePars,fluorophorePosPx,fluorophorePosPxYOYO);
end

figure,imagesc(finalImageSCMOS(:,:,1))
figure,imagesc(finalImageBGYOYO)

import Simulate.save_image_separate;
%  [name1,name2] = save_image_separate(noisyImage2D,noisyImageYoyo,nameext)
[name1,name2] = save_image_separate(finalImageSCMOS,finalImageBGYOYO,'');


% TEST.sdd_run_tests('output')
