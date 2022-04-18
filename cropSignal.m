function croppedSignal = cropSignal(signal, windowLength, numOfAddedSamples)

% Crop left and right extensions and return the right portion of signal

% 2021, Sepehr Ghanbari, Pavel Rajmic

croppedSignal = signal(windowLength+1:end-windowLength-numOfAddedSamples);

% % Sepehr's orig code:
% %% Removing Extra Zeros
% extraWinsZeros = 2 * progressStep;
% signalAppendedZeros = length(finalSignal) - extendedSignalLen - 2 * extraWinsZeros;
% originalBegin = ceil(signalAppendedZeros/2) + extraWinsZeros + 1;
% originalEnd = originalBegin + extendedSignalLen - 1;
% finalSignal = finalSignal(originalBegin:originalEnd,1);
% 