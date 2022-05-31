%function to minimize on W and Z
%it is necessary to evaluate it repeatedly for obtaining step sizes
function value = fGammaEvaluate(W,Z,gamma)
    value = gamma*(norm(W,'fro')^2) + (norm(transpose(Z)*W,'fro')^2);
end

%W=sqrt(epsilon)*fZ
%Z=transpose(fZ)

% epsilon*transpose(fZ)*sqrt(epsilon)*fZ