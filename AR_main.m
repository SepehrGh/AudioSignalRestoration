clear all
close all
clc

% requires LTFAT to generate the Hann window
%% Input parameters
% a = [0.7 0.25 -0.8];  %AR coefficients
% d = [0.3 0.9 1];      %boundary conditions (same length as the coefs)
% inputSignalLength = 1500;   %length of signal
% v = 0.01;            %variance of noise (no noise <=> v=0)

%% Input signal (generate AR process)
% [originalSignal,~] = myAR(a, d, v, inputSignalLength);
% originalSignal = originalSignal/abs(max(originalSignal)); %Not sure if normalization has to be done?


% put the clipped audio samples in this path (SurveySamples\Clipped\)
% and the original ones here (SurveySamples\Original\)
% after running the code the results will be in the current directory under "Declipped\"

resultsDir = "Declipped\";
mkdir(resultsDir)

filenameList = ["a08_violin"];
%"a42_accordion", "bassoon", "a18_bassoon", "a41_celesta", "a16_clarinet", "a35_glockenspiel", 
%"a58_guitar_sarasate", "a25_harp", "a60_piano_schubert", "a08_violin", "a66_wind_ensemble_stravinsky"
for i=1:length(filenameList)
    filename = filenameList(i);
    %% Read original signal (unclipped)
    originalDir = "SurveySamples\Original\" ... 
        + filename + ".wav";
    [originalSignal,~] = audioread(originalDir);
    %% Loop on inputSDR values
    inputSdrList = ["05"]; %"01dB", "03dB", "07dB", "10dB", "15dB", "20dB"
    for j=1:length(inputSdrList)
        inputSdr = inputSdrList(j);
        finalPath = resultsDir + inputSdr + "dB\";
        mkdir(finalPath)
        %% Read clipped signal
        clippedDir = "SurveySamples\Clipped\" ...
           + inputSdr + "dB\" + filename + "_" + inputSdr + ".wav";
        [clippedSignal,Fs] = audioread(clippedDir);
        inputSignalLength = length(clippedSignal);
%         clippingThreshold = .6; %0.98

        %% Extract Clipping Level
%         [clippedSignal, clippedIndexesLow, clippedIndexesHigh]  = ...
%             clipSignal(originalSignal, -clippingThreshold, clippingThreshold);
        [lowThreshold, highThreshold, clippedIndexesLow, clippedIndexesHigh] = ...
            extractClippingLevel(clippedSignal);
        clippingThreshold = max(abs(lowThreshold), highThreshold);
        initialClippedIndexes = or(clippedIndexesLow, clippedIndexesHigh); %union of indexes

        %% Plot
        figure;
%         time-domain signals
        plot([originalSignal clippedSignal]);
        title('Input signal and its clipped version');
        xlabel('time (samples)');

        %and mark clipped samples
    %     hold on
    %     clippedIndexesLowJustThem = 1:inputSignalLength;
    %     clippedIndexesLowJustThem = clippedIndexesLowJustThem(clippedIndexesLow);
    %     stem(clippedIndexesLowJustThem,zeros(length(clippedIndexesLowJustThem),1),'g')
    %     clippedIndexesHighJustThem = 1:inputSignalLength;
    %     clippedIndexesHighJustThem = clippedIndexesHighJustThem(clippedIndexesHigh);
    %     stem(clippedIndexesHighJustThem,zeros(length(clippedIndexesHighJustThem),1),'m')

        %% Levinson Estimation
        %[levinsonCoeffs,~] = generateLevinsonEstimation(a,s);
        %disp(levinsonCoeffs);
        %fprintf('Norm difference with Levinson esimation: %f\n',norm(levinsonCoeffs-a));

        %% Hankel Matrix
        %H = generateHankelMatrix(x, 121, 190);
        %fprintf('Hankel matrix rank: %j\n\n', rank(H));
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

        %split into chunks, window
        [chunks, chunkIndices, chunksClippedIndicesLow, chunksClippedIndicesHigh, analysisWindow] = ...
            applyAnalysisWindow(windowType, windowLength, hopSize, extClippedSignal, clippedIndexesLow, clippedIndexesHigh);
        totalNumOfWindows = size(chunkIndices,2);

        %% Processing chunks using NSAO
        lowRankChunks = zeros(windowLength,totalNumOfWindows);  %allocate space
        fGammaValuesRecords = cell(totalNumOfWindows,1);
        chunksClippedIndicesLow = logical(chunksClippedIndicesLow);
        chunksClippedIndicesHigh = logical(chunksClippedIndicesHigh);
        clippedIndexes = or(clippedIndexesLow, clippedIndexesHigh); %union of indexes
        waithandle = waitbar(0,'calculating now');
        for s = 1:100 %totalNumOfWindows
            if (all(chunksClippedIndicesLow(:,s) == 0) && all(chunksClippedIndicesHigh(:,s) == 0)) %if there is no clipped sample in chunk, skip chunk entirely
                lowRankChunks(:,s) = chunks(:,s);
                continue
            end
            Z = generateHankelMatrix(chunks(:,s));
            [lowRankZ, fGammaValues] = ...
                nsaoGpm(Z, 0.01, 1, chunksClippedIndicesLow(:,s), chunksClippedIndicesHigh(:,s), clippingThreshold);
            %gamma = 10^(-2), eta = 1 (previously 1.1)
            firstCol = lowRankZ(:,1);
            lastRow = lowRankZ(end,2:end);
            lowRankChunks(:,s) = vertcat(firstCol,transpose(lastRow));
            fGammaValuesRecords{s} = fGammaValues;  %store the current values
            progress = s/totalNumOfWindows;
            waitbar(progress,waithandle);
            disp(progress*100);
        end
        close(waithandle)
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

        %% SDRc
        deltaSDRc = sdrC(originalSignal, finalSignal, initialClippedIndexes) - sdrC(originalSignal, clippedSignal, initialClippedIndexes)

        %% Write result
    %     dir = "C:\Users\Sepehr\Desktop\Audio Restoration\newImpl (21.3.22)\declipped\05dB\" + filename + "\";
    %     writeFilename = dir + "clipped_CL_" + clippingThreshold + ".wav";
    %     audiowrite(writeFilename, clippedSignal, Fs);
    %     writeFilename = dir + "CL_" + clippingThreshold + "_HS_" + windowLength +  ".wav";
    %     audiowrite(writeFilename, finalSignal, Fs);

        hankelSize = (windowLength+1)/2;
        writeFilename = finalPath + filename + "_HankelSize_" + hankelSize +  ".wav";
        audiowrite(writeFilename, finalSignal, Fs);
        %% Plots
        figure

        plot([originalSignal clippedSignal finalSignal])
        % hold on
        xlabel('time (samples)')
        title('Original, clipped and recovered signals')
        legend('original','clipped','recovered')
        axis tight

%         plot([originalSignal finalSignal])
%         % hold on
%         xlabel('time (samples)')
%         title('Original and recovered signals')
%         legend('original', 'recovered')
%         axis tight


%         %% Plot Seperately
%         figure
%         subplot(2,1,1)
%         plot(clippedSignal)
%         xlabel('time (samples)')
%         % ylabel('Amplitude')
%         title('Input Signal')
%         
%         subplot(2,1,2)
%         plot(finalSignal);
%         xlabel('time (samples)')
%         % ylabel('Amplitude')
%         title('Final Signal')

        figure
        plot([originalSignal-finalSignal])
        xlabel('time (samples)')
        title('Difference between original and recovered signals')
        axis tight

        figure
        for s = 1:totalNumOfWindows
            plot(fGammaValuesRecords{s})
            hold on
        end
        title('Value of the objective function during half-iterations')
    end
end