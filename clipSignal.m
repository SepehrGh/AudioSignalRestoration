function [clippedSignal, clippedIndexesLow, clippedIndexesHigh] = ...
    clipSignal(signal, lowThreshold, highThreshold)

% Clip the signal using thresholds

% 2021, Sepehr Ghanbari, Pavel Rajmic

clippedSignal = signal;

clippedIndexesLow  = (signal < lowThreshold);
clippedIndexesHigh = (signal > highThreshold);
clippedSignal(clippedIndexesLow)  = lowThreshold;
clippedSignal(clippedIndexesHigh) = highThreshold;

end