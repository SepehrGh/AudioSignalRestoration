clear all
close all
clc

% requires LTFAT to generate the Hann window

%% Input parameters
a = [0.7 0.25 -0.8];  %AR coefficients
d = [0.3 0.9 1];      %boundary conditions (same length as the coefs)
inputSignalLength = 1500;   %length of signal
v = 0.01;            %variance of noise (no noise <=> v=0)

clippingThreshold = .6; %0.98

%% Input signal (generate AR process)
[inputSignal,~] = myAR(a, d, v, inputSignalLength);
inputSignal = inputSignal/abs(max(inputSignal)); %Not sure if normalization has to be done?

% inputSignal = load('inputSignal1.mat');
% inputSignal = [inputSignal.inputSignal];
%Running with sinusoid signal as input
% t=1:1:500;
% f=10;
% x=sin((1/500)*2*pi*f*t);
% inputSignal = transpose(x);

%% Signal Clipping
[clippedSignal, clippedIndexesLow, clippedIndexesHigh]  = ...
    clipSignal(inputSignal, -clippingThreshold, clippingThreshold);
initialClippedIndexes = or(clippedIndexesLow, clippedIndexesHigh); %union of indexes

%% Plot
figure;
%time-domain signals
plot([inputSignal clippedSignal]);
title('Input signal and its clipped version');
xlabel('time (samples)');
%and mark clipped samples
hold on
clippedIndexesLowJustThem = 1:inputSignalLength;
clippedIndexesLowJustThem = clippedIndexesLowJustThem(clippedIndexesLow);
stem(clippedIndexesLowJustThem,zeros(length(clippedIndexesLowJustThem),1),'g')
clippedIndexesHighJustThem = 1:inputSignalLength;
clippedIndexesHighJustThem = clippedIndexesHighJustThem(clippedIndexesHigh);
stem(clippedIndexesHighJustThem,zeros(length(clippedIndexesHighJustThem),1),'m')

%% Levinson Estimation
%[levinsonCoeffs,~] = generateLevinsonEstimation(a,s);
%disp(levinsonCoeffs);
%fprintf('Norm difference with Levinson esimation: %f\n',norm(levinsonCoeffs-a));

%% Hankel Matrix
%H = generateHankelMatrix(x, 121, 190);
%fprintf('Hankel matrix rank: %i\n\n', rank(H));
%[~,S,~] = svd(H);

%% Applying Analysis Window
windowType = 'rect';  %window type
windowLength = 99;  %only odd lengths allowed; otherwise Hankel matrices could not be formed
% overlap = 22;
hopSize = 33;   %time step (hop size) of the window

%extend signal to make it suitable for analysis
[extClippedSignal, numOfAddedSamples] = ...
    extendSignal(clippedSignal, windowLength, hopSize);
%extend also the masks
[clippedIndexesLow, ~] = ...
    extendSignal(clippedIndexesLow, windowLength, hopSize);
clippedIndexesLow = logical(clippedIndexesLow);
[clippedIndexesHigh, ~] = ...
    extendSignal(clippedIndexesHigh, windowLength, hopSize);
clippedIndexesHigh = logical(clippedIndexesHigh);
clippedIndexes = or(clippedIndexesLow, clippedIndexesHigh); %union of indexes

%split into chunks, window
[chunks, chunkIndices, chunksClippedIndicesLow, chunksClippedIndicesHigh, analysisWindow] = ...
    applyAnalysisWindow(windowType, windowLength, hopSize, extClippedSignal, clippedIndexesLow, clippedIndexesHigh);
totalNumOfWindows = size(chunkIndices,2);

%% Processing chunks using NSAO
lowRankChunks = zeros(windowLength,totalNumOfWindows);  %allocate space
fGammaValuesRecords = cell(totalNumOfWindows,1);
% start = 1;
% len = (windowLength+1)/2;
chunksClippedIndicesLow = logical(chunksClippedIndicesLow);
chunksClippedIndicesHigh = logical(chunksClippedIndicesHigh);
for s = 1:totalNumOfWindows
    if (all(chunksClippedIndicesLow(:,s) == 0) && all(chunksClippedIndicesHigh(:,s) == 0)) %if there is no clipped sample in chunk, skip chunk entirely
        lowRankChunks(:,s) = chunks(:,s);
        continue
    end
    Z = generateHankelMatrix(chunks(:,s));
    [lowRankZ, fGammaValues] = ...
        nsaoGpm(Z, 0.01, 1.1, chunksClippedIndicesLow(:,s), chunksClippedIndicesHigh(:,s), clippingThreshold); %gamma = 10^(-2) (older version:0.8), eta = 1.1
    firstCol = lowRankZ(:,1);
    lastRow = lowRankZ(end,2:end);
    lowRankChunks(:,s) = vertcat(firstCol,transpose(lastRow));
    fGammaValuesRecords{s} = fGammaValues;  %store the current values
%     fprintf('Matrix rank: %i\n\n', rank(Z));
%     fprintf('Matrix rank after NSAO-GPM: %i\n', rank(lowRankZ));
end

%% Applying Synthesis Window
windowType = 'hann';
% windowLength = 33;
% overlap = 22;
% hopSize = 11;
[synthSignal,synthesisWindow] = ...
    applySynthesisWindow(windowType, hopSize, lowRankChunks, chunkIndices, length(extClippedSignal), analysisWindow);

%Cropping the right and left artificial extensions
finalSignal = cropSignal(synthSignal, windowLength, numOfAddedSamples);

% %% Scaling Final Signal
% scalingFactor = max(clippedSignal) / max(synthSignal);
% synthSignal = synthSignal .* scalingFactor;

%% Evaluate numerically
norm(finalSignal-inputSignal);

%% SDRc
deltaSDRc = sdrC(inputSignal, finalSignal, initialClippedIndexes) - sdrC(inputSignal, clippedSignal, initialClippedIndexes)

%% Plots
figure
plot([inputSignal clippedSignal finalSignal])
% hold on
% plot(clippedSignal)
xlabel('time (samples)')
% ylabel('Amplitude')
title('Original, clipped and recovered signals')
legend('original','clipped','recovered')
axis tight

% %% Plot Seperately
% figure
% subplot(2,1,1)
% plot(clippedSignal)
% xlabel('time (samples)')
% % ylabel('Amplitude')
% title('Input Signal')
% 
% subplot(2,1,2)
% plot(finalSignal);
% xlabel('time (samples)')
% % ylabel('Amplitude')
% title('Final Signal')

figure
plot([inputSignal-finalSignal])
xlabel('time (samples)')
% ylabel('Amplitude')
title('Difference between original and recovered signals')
axis tight

figure
for s = 1:totalNumOfWindows
    plot(fGammaValuesRecords{s})
    hold on
end
title('Value of the objective function during half-iterations')
