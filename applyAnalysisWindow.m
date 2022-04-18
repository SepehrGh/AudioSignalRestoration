function [chunks,indices,chunksClippedIndicesLow,chunksClippedIndicesHigh,window] = ...
    applyAnalysisWindow(type, winLen, hopSize, signal, clippedIndexesLow, clippedIndexesHigh)

% Splits the signal into possibly overlapping segments (chunks) and weights
% them according to the window type

% 2021, Sepehr Ghanbari, Pavel Rajmic

%% Check the window type
if type == 'rect'
    window = ones(winLen,1);
elseif type == 'hann'
    window = fftshift(firwin(g,w)); %not tested; uses LTFAT toolbox
    %!!! should any window other than rect be used, then also the
    %clipping threshold should be windowed and passed further
else
    error('Unknown type of window')
end

%% Calculating number of windows
signalLen = length(signal);
% totalNumOfWindows = ((length(signal) - winLen) / hopSize);
totalNumOfWins = signalLen / hopSize;

%% Chunking with overlap
indices = zeros(winLen,totalNumOfWins);   %preallocate matrix where indices of windows will be stored
chunks  = zeros(winLen,totalNumOfWins);   %preallocate matrix where windowed signal segments will be stored
chunksClippedIndicesLow  = zeros(winLen,totalNumOfWins); %preallocate matrix where clipping indices will be stored
chunksClippedIndicesHigh = zeros(winLen,totalNumOfWins); %preallocate matrix where clipping indices will be stored

%indices of windows (will be used in synthesis as well)
firstWindowIndices = 1 - floor(winLen/2) : ceil(winLen/2); %first window, nonperiodized
for s = 1:totalNumOfWins
    indices(:,s) = firstWindowIndices + (s-1)*hopSize;
end
indices = 1 + mod(indices-1, signalLen);  %modulo => periodized system

% In the following loop, input signal is divided into chunks of size
% "winLen", with regards to the overlap specified, and windowed
for s = 1:totalNumOfWins
    chunks(:,s) = window .* signal(indices(:,s));  %window the signal chunks
    chunksClippedIndicesLow(:,s)  = clippedIndexesLow(indices(:,s)); %book the clipped indices within chunks
    chunksClippedIndicesHigh(:,s) = clippedIndexesHigh(indices(:,s)); %book the clipped indices within chunks
end

% %% Preparing Weighing Window
% if type == 'hann'
%     window = hann(winLen);
% end
% %TODO: add other possible window types
% %% Normalizing Window
% window = normalize(window,'norm',1); %Empirically proven for Hann Window
%
% %     General way of normalizing the window
% %     sum = zeros(1,length);
% %     for i = 1:length
% %         sum = sum + window(i);
% %     end
% %     window = window./transpose(sum);
%
% %% Test Partition of Unity
% sum = zeros(1,winLen);
% for i = 1:winLen
%     sum = sum + window(i);
% end
% %% Weighing Chunks
% for i = 1:totalNumOfWins
%     chunks{i} = chunks{i}.*window;
% end

