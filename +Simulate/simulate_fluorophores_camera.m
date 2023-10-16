function [finalImage,final1D,fluorophorePosPx,finalImageBG,final1DBG,finalImageBGYOYO,finalImageBGYOYO1D] =...
    simulate_fluorophores_camera(camera,sim,chipPars,labelPars,noisePars,fluorophorePosPx, fluorophorePosPxYOYO)
    % function to simulate fluorophores in 2D
    % 
    %     simulate_fluorophores_camera
    

    % fluorophore positions
    numberPeaks = sim.numberPeaks;
    gapBetweenPeaks = sim.gapBetweenPeaks;
    distTwinPeaks = sim.distTwinPeaks;
    stgap = sim.stgap;
    nY = sim.nY;
    nX = sim.nX;
    
    if nargin < 6
        % fixed evenly distributed fluorophore positions
        [fluorophorePos, fluorophorePosPx] = fixed_fluorophores(numberPeaks, gapBetweenPeaks, distTwinPeaks,stgap, labelPars.nmbp, chipPars.pixelsize, nY);
    end
    
    % photon hitting probabilities / for single fluorophorePosPx
    [pxPhotons, img] = gen_photon_prob(noisePars.lambdaSg, noisePars.lambdaBg, fluorophorePosPx, noisePars.sigma, nX, nY);

    
    % here we consider the cases for scmos and emccd.
    switch camera
        case 'scmos'
            noisyImage = noise_model_scmos(pxPhotons', chipPars);
        case 'emccd'
            % noise model EMCCD
            noisyImage = noise_model_emccd(pxPhotons', chipPars);
        otherwise
            
    end
    
    % final image
    finalImage = reshape(noisyImage,nY,nX);
    % final image 1D - only center row
    final1D = finalImage(round(nY/2),:);
    
    if nargout >=4
        photonsPlacedBG = noisePars.lambdaBg*ones(length(pxPhotons),1);
        switch camera
            case 'scmos'
                noisyImageBG = noise_model_scmos(photonsPlacedBG', chipPars);
            case 'emccd'
                % noise model EMCCD
                noisyImageBG = noise_model_emccd(photonsPlacedBG', chipPars);
            otherwise    
        end
        finalImageBG = reshape(noisyImageBG,nY,nX);
        final1DBG = finalImageBG(round(nY/2),:);
    end
    
    if nargout >=6
        % also yoyo (if narguot asks for it only)
        if nargin < 7
            try % place yoyo every pixel (so roughly 1 yoyo per 500 bp
            fluorophorePosPxYOYO = [[(fluorophorePosPx(1,1)-15):1:(fluorophorePosPx(end,1)+15)]' round(nY/2)*ones(length(fluorophorePosPx(1,1)-15:1:fluorophorePosPx(end,1)'+15),1)];
            catch
            fluorophorePosPxYOYO = [[(fluorophorePosPx(1,1)-5):1:(fluorophorePosPx(end,1)+5)]' round(nY/2)*ones(length(fluorophorePosPx(1,1)-5:1:fluorophorePosPx(end,1)'+5),1)];
            end % since there should be roughly 1 yoyo per 4 bp, we have 100 times lambdaS. lambdaYoyo
        end
        [pxPhotonsYOYO] = gen_photon_prob(100*noisePars.lambdaYoyo,noisePars.lambdaBg,fluorophorePosPxYOYO,noisePars.sigma, nX, nY);
%         pxPhotonsYOYO = pxPhotonsYOYO/sum(pxPhotonsYOYO(:));
%         [photonsPlacedYOYO]= place_photons(pxPhotonsYOYO,noisePars.lambdaSg*length(fluorophorePos)+noisePars.lambdaBg*length(pxPhotons));
        switch camera
            case 'scmos'
                noisyImageYOYO = noise_model_scmos(pxPhotonsYOYO', chipPars);
            case 'emccd'
                % noise model EMCCD
                noisyImageYOYO = noise_model_emccd(pxPhotonsYOYO', chipPars);
            otherwise    
        end

        finalImageBGYOYO = reshape(noisyImageYOYO,nY,nX);
        finalImageBGYOYO1D = finalImageBGYOYO(round(nY/2),:);

    end
    

end

function [fluorophorePos, fluorophorePosPx] = fixed_fluorophores(npeaks,ngapPeaks,dpeaks,stgap,nmbp,pixelsize,nY)
    % npeaks - number of peaks
    % ngappeaks -difference between leftmost peaks
    % dpeaks - difference to second peak
    % here we get random fluorophore positions
%     npeaks = 20; % number of peaks/fixed
%     ngapPeaks = 15000; % gap between leftmost peaks /fixed
    if nargin < 4
        stgap = 20000;
    end
    
    firstPos = stgap+ngapPeaks:ngapPeaks:(stgap+npeaks*ngapPeaks); 
    secondPos=  firstPos+dpeaks; % dpeaks tunable
    extra = [];
%     extra = secondPos(end)+[1000 2000 4000 8000 16000];
    fluorophorePos = sort([firstPos secondPos extra])';
%     rmap = zeros(firstPos(end)+ngapPeaks,1);
%     rmap(firstPos) = 1;
%     rmap(secondPos) = 1;
%     rmap2 =  zeros(firstPos(end)+ngapPeaks,1);
%     rmap2(firstPos(1)-2000:1000:secondPos(end)+2000)=1;
    fluorophorePosPx = [fluorophorePos*nmbp/pixelsize round(nY/2)*ones(length(fluorophorePos),1)];

end

function [pxPhotons,img] = gen_photon_prob(lambdaSig,lambdaBg,fluorophorePosPx,sigma,nX,nY)
% lambdaSig = 3000;
% lambdaBg = 100;
% sigma =1.5;

% based on
% Huang F, Schwartz SL, Byars JM, Lidke KA. Simultaneous multiple-emitter fitting for single molecule super-resolution imaging. Biomedical optics express. 2011 May 1;2(5):1377-93.

% beadPos = [50.2 40.3];

% nX = 100;
% nY = 100;
[X,Y] = meshgrid(1:nX,1:nY);

X = X(:);
Y = Y(:);
pxPhotons = zeros(size(X));

deltaEx = @(x,b) 1/2*(erf((x-b+1/2)/(sqrt(2)*sigma))-erf((x-b-1/2)/(sqrt(2)*sigma)));

for i=1:size(fluorophorePosPx,1)
    pxPhotons = pxPhotons+lambdaSig*deltaEx(X,fluorophorePosPx(i,1)).*deltaEx(Y,fluorophorePosPx(i,2));
end
pxPhotons = pxPhotons+lambdaBg;
% pxPhotons = pxPhotons/sum(pxPhotons(:));
img = reshape(pxPhotons,nY,nX);


end
% 
% function [imgRow]= place_photons(probs,nPhotons)
% 
% 
%     nPixels = length(probs); % number of pixels
% %     nPhotons = 30;
%     % Choose to distribute the photons over the pixels
%     %% DO we need this or we directly calculate nic, etc?
%     randPixels = randsrc(nPhotons,1,[1:nPixels;probs']);
%     [cnt_unique, unique_a] = hist(randPixels,unique(randPixels));
% 
% 
%     imgRow = zeros(1,nPixels);
%     imgRow(unique_a) = cnt_unique; % here 
% 
% %     photMat = reshape(imgRow,length(gridX),length(gridY));
% 
% end

% noise image emccd
function noisyImage=noise_model_emccd(photonsPlaced,chipPars)
    %
    % from Super-resolution fight club: assessment  of 2D and 3D single-molecule 
    %localization microscopy software
    gain = chipPars.gain;
    roNoise = chipPars.roNoise;
    adFactor = chipPars.adFactor;
    offset = chipPars.offset;
    c = chipPars.c; 
    QE = chipPars.QE;
 

    poissVals = poissrnd(photonsPlaced*QE+c);
    % Generate amplified electric signal
    elecImage = gamrnd(poissVals,gain);

    % Perform readout on each pixel     
    noisyImage = (roNoise/adFactor) * randn(1,length(photonsPlaced)) + offset + elecImage/adFactor;
    
    % round (uniformly distributed noise)
    noisyImage = round(noisyImage);

end

% noise image scmos
function noisyImage=noise_model_scmos(photonsPlaced, chipPars)
    %
    % Physics-Based Noise Modeling for Extreme
    % Low-Light Photography
    gain = chipPars.gain;
    roNoise = chipPars.roNoise;
%     adFactor = chipPars.adFactor;
    offset = chipPars.offset;
%     c = chipPars.c; 
%     QE = chipPars.QE;
 

    poissVals = gain*poissrnd(photonsPlaced);
    % Generate amplified electric signal
%     elecImage = gamrnd(poissVals,gain);
% random('logistic', 0, sigma_Re) 
    Nread = random('logistic', 0, roNoise,1,length(photonsPlaced));
%     Nread = rand(1,length(photonsPlaced));
%     Nread = roNoise.*log(Nread./(1-Nread));

    % Perform readout on each pixel     random('logistic', 0, sigma_Re)
    noisyImage = poissVals + offset + Nread;
    % round (uniformly distributed noise)
    noisyImage = round(noisyImage);

end
