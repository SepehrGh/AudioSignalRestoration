function [lowRankZ, fGammaValues] = nsaoGpm(Z, gamma, eta, clippedIndexesLow, clippedIndexesHigh, clippingThreshold) 

%implementation of the NSAO-GPM algorithm
%gamma>0, eta>=1

% 2021, Sepehr Ghanbari, Pavel Rajmic

dimensionsZ = size(Z); %[mZ,nZ] (mZ <= nZ)
wVector = ones(1,min(dimensionsZ)); %vector of ones of length mZ
W = diag(wVector); %initial matrix W - mZ*mZ (main diagonal = 1)

clippedIndexes = or(clippedIndexesLow, clippedIndexesHigh); %union of indexes
% gradientMask = clippedIndexes; % Mask in order to have zeroes on the unclipped samples in the gradient
gradientMask = generateHankelMatrix(clippedIndexes);

fGamma = fGammaEvaluate(W,Z,gamma); %function to minimize w.r.t. W and Z
fGammaValues = [fGamma];

epsilon = 1e-8;
counter = 0;
maxIterations = 100;
while 1
    if all(Z(:) == 0)
        break
    end
    if counter > maxIterations
        break
    end
%     if length(fGammaValues)>1
%         currfGamma = fGammaValues(end);
%         prevfGamma = fGammaValues(end-1);
% %         if (abs((currfGamma - prevfGamma)/prevfGamma)<(10^(-5))) %if the amount of decrease is tiny, finish running algorithm
%         if (norm((currfGamma - prevfGamma),'fro')<(10^(-5))) %if the amount of decrease is tiny, finish running algorithm
%             disp(counter)
%             break
%         end
%     end

%     dW = 2*(gamma*W + Z*transpose(Z)*W);
    dW = (gamma*W + Z*transpose(Z)*W);
    fW = calculateProjectionW(dW) - eye();
    alphaW = ((-1)*trace(transpose(fW)*dW))/(2*fGammaEvaluate(fW,Z,gamma));
%     alphaW = ((-1)*norm(fW,'fro')^2)/fGammaEvaluate(fW,Z,gamma);
    W = W + alphaW*fW;
    
    fGammaValues = [fGammaValues fGammaEvaluate(W,Z,gamma)];  %record the value
    
    dZ = transpose(Z)*(W*transpose(W) + gamma*epsilon*eye());
    fZ = dZ .* gradientMask;
    fZ = calculateProjectionZ(fZ);
    alphaZ = ((-1)*trace(transpose(fZ)*dZ))/(2*(norm(transpose(fZ)*W,'fro')^2)); 
    Z = Z + alphaZ*fZ;
    Z = applyClippingConstraint(Z, clippedIndexesLow, clippedIndexesHigh, clippingThreshold);
    gamma = gamma / eta; %gamma = max(gamma/eta,gamma_min);
    
    fGammaValues = [fGammaValues fGammaEvaluate(W,Z,gamma)];  %record the value
    
    counter = counter + 1;
end
lowRankZ = Z;