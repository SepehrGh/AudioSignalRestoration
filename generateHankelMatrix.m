function hankelMatrix = generateHankelMatrix(inputSignal) %, firstColumnStarting, firstColumnFinishing, lastRowFinishing)

%generates a square Hankel matrix from a time-domain vector

% 2021, Sepehr Ghanbari, Pavel Rajmic

signalLen = length(inputSignal);
if rem(signalLen,2) == 0  %is even?
    error('Even lengths not allowed');
end

switchIndex = (signalLen+1)/2;
hankelMatrix = hankel(inputSignal(1:switchIndex),inputSignal(switchIndex:signalLen));

% firstColumn = inputSignal(firstColumnStarting:firstColumnFinishing);
% lastRow = inputSignal(firstColumnFinishing:lastRowFinishing);
% hankelMatrix = hankel(firstColumn,lastRow);