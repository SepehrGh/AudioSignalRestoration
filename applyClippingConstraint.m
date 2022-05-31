function [res] = applyClippingConstraint(mat, clippedIndexesLow, clippedIndexesHigh, clippingThreshold)
%     clippedIndexesLow = logical(clippedIndexesLow);
%     clippedIndexesHigh = logical(clippedIndexesHigh);
    firstCol = mat(:,1);
    lastRow = mat(end,2:end);
    vec = vertcat(firstCol,transpose(lastRow));
    res = vec;
    indexesToLower = and(clippedIndexesLow, vec > (-clippingThreshold));
    indexesToRaise = and(clippedIndexesHigh, vec<clippingThreshold);
    res(indexesToLower) = -clippingThreshold;
    res(indexesToRaise) = clippingThreshold;
    
%% Generating Output Hankel Matrix
    res = generateHankelMatrix(res);
end