function [extendedSignal, numOfAddedSamples] = extendSignal(signal,winLen,hopSize)

% 2021, Sepehr Ghanbari, Pavel Rajmic

% A/ to make the system periodic (signalLen divisible by hopSize)
signalLen = length(signal);
numOfAddedSamples = hopSize*ceil(signalLen/hopSize) - signalLen; %how many samples should be appended...
extendedSignal = [signal; zeros(numOfAddedSamples,1)];                   %...to make the length divisible by hopSize

% B/ to avoid effects from periodized windows
extendedSignal = [zeros(winLen,1); extendedSignal; zeros(winLen,1)]; %extend; still divisible by a; one could find shorter appendices but this is simple
% signalLen = length(extendedSignal);

% Sepehr's original code:
% %% Appending initial and final zeros
% %Adding 2*progressStep length of zeros at the beginning and the end
% %in order to have the required number of windows at the beginning and
% %also at the end of the signal.
% numOfAppendedWindows = ceil((winLen / hopSize) - 1);
% numOfZeros = numOfAppendedWindows * hopSize;
% signal = vertcat(zeros(numOfZeros+ceil(numOfAddedSamples/2),1), signal);
% signal = vertcat(signal, zeros(numOfZeros+floor(numOfAddedSamples/2),1));
