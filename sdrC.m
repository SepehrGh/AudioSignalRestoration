function [sdrC] = sdrC(sig1,sig2,clippedIndexes)

%extract clipped section
sig1C = sig1(clippedIndexes);
sig2C = sig2(clippedIndexes);

%calculate sdrC
sdrC = 20*log(norm(sig1C)/norm(sig1C-sig2C));

end