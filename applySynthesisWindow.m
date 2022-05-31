function [synthSignal, synthesisWindow] = ...
    applySynthesisWindow(winType, hopSize, chunks, indices, extendedSignalLen, analysisWindow)

% Synthesizes the possibly overlapping signal segments (chunks) using the synthesis window

% 2021, Sepehr Ghanbari, Pavel Rajmic

%% Preparing Synthesis Window
winLen = size(indices,1);
if winType == 'hann'
%     window = hann(winLen);
      synthesisWindow = fftshift(firwin('hann',winLen));
end
%check whether the windows are dual and get the scaling factor
[scal,residual] = gabdualnorm(analysisWindow,synthesisWindow,hopSize,winLen);
scal = real(scal); %Quick fix
if residual > 10e-10
    error('Analysis and Synthesis windows are probably not dual!')
end
synthesisWindow = synthesisWindow*winLen/scal;  %scaling to obtain proper amplitude at synthesis

% window = normalize(window,'norm',1); %Empirically proven for Hann Window
%     General way of normalizing the window
%     sum = zeros(1,length);
%     for i = 1:length
%         sum = sum + window(i);
%     end
%     window = window./transpose(sum);

%% Reconstructing Final Signal Using Chunks
% progressStep = winLen - overlap;
totalNumOfWindows = size(indices,2);
synthSegments = zeros(size(indices));    %preallocate for the synthesized chunks
synthSignal = zeros(extendedSignalLen,1); %preallocate space for the synthesized signal
% finalSignal = zeros(progressStep*(4+(totalNumOfWindows-1)),1); %extra zeros: 2 progressWindows at the beginning and 2 at the end
%     numOfAppendedWindows = ceil((winLen / progressStep) - 1);
%     numOfZeros = numOfAppendedWindows * progressStep; %appended both at the beginning and the end
%     finalSignal = zeros(signalLen+2*numOfZeros,1);
%     totalNumOfWindows = numOfWindows + 2*numOfAppendedWindows;
for s=1:totalNumOfWindows
    synthSegments(:,s) = chunks(:,s) .* synthesisWindow;  %apply synthesis window to each chunk
    synthSignal(indices(:,s)) = synthSignal(indices(:,s)) + synthSegments(:,s);  %add the current synthesis
%     winProd = signal{s}.*synthesisWindow;
%     beginIndex = (s-1)*progressStep + 1;
%     endIndex = beginIndex + winLen - 1;
%     finalSignal(beginIndex:endIndex) = finalSignal(beginIndex:endIndex) + winProd;
end
