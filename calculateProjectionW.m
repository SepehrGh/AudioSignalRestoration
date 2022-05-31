function [res] = calculateProjectionW(matrix)
    res = matrix;
    d = logical(eye(size(matrix)));
    res(d) = 1;  %substitute the diagonal with ones
end
