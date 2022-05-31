% function [res] = calculateProjectionZ(matrix, clippedIndexesLow, clippedIndexesHigh, clippingThreshold) %matrix = dZ
function [res] = calculateProjectionZ(matrix) %matrix = dZ
%projects the matrix to the space of Hankel matrices and then on the set of consistent solutions
%(note that this is most probably not the same as projecting onto
% the intersection of the sets, but should suffice in practice)
%
%matrix can be rectangular (BUT CURRENT CODE BY PAVEL ONLY WORKS FOR SQUARE MATRICES)

% 2021, Sepehr Ghanbari, Pavel Rajmic

% matrix = generateHankelMatrix(1:33) + randn(size(matrix));  %just for testing purposes

% %% Calculating fs Values
% r = size(matrix,1)-1; %M_rows = r+1 (0,1,...,r)
% N = (size(matrix,2)+r)-1; %M_cols = N-r+1 (0,1,...,N-r) -> N = M_cols+r-1
% c = N-r;
% f = zeros(N+1,1); %f0,f1,...,fN
% for s = 0:N %0,1,2,...,N
%     iS = [];
%     for i = 0:r
%         for j = 0:c %c = N-r (# of columns in dZ) (To be considered: InnerLoopCount=1350 for 10*10 matrix)
%             if i+j == s
%                 iS = [iS matrix(i+1,j+1)];
%                 break
%             end
%         end
%     end
%     f(s+1) = (1/length(iS))*sum(iS);  %taking the mean
%     %% Applying Clipping Constraint
%     %         f(s+1) = applyClippingConstraint(s+1, f(s+1), clippedIndexes, clippingThreshold);
% end

extractedDiagonals = fliplr(spdiags(fliplr(matrix),-size(matrix,1)+1:size(matrix,2)-1)); %extracts all antidiagonals from the Hankel matrix
columnSum = sum(extractedDiagonals); %columnwise summation
divisors = [1:size(extractedDiagonals,1) size(extractedDiagonals,1)-1:-1:1];  %this should be generalized for the rectangular case
f = (columnSum ./ divisors)'; %divide, i.e. get the columnwise means
% norm(f-ff)  %should be very close to zero

%% Applying Clipping Constraint
% f = applyClippingConstraint(f, clippedIndexesLow, clippedIndexesHigh, clippingThreshold);

%% Generating Output Hankel Matrix
res = generateHankelMatrix(f); %,1,r+1,N+1);

