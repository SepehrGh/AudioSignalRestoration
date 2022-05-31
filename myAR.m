function [x,noise] = myAR(a, x0, v, N)

% Generate AR process using the definition (i.e., recursively)
%
% x ... AR process
% noise ... driving noise that was used for the generation
%
% a ... vector of AR coefs
% x0 ... vector of boundary conditions
% N ... length of the process (output signal)
% v ... variance of the Gaussian noise (sigma squared)

% 2021, Sepehr Ghanbari, Pavel Rajmic

p = length(x0);               % learn the model order from boundary conditions
noise = sqrt(v)*randn(N+p,1); % generate Gaussian driving noise (first p entries are useless)

x = zeros(N+p,1);            % allocate vector
x(1:p) = x0;                 % copy boundary conditions as the beginning of the output

for j = 1:N                % loop over output signal elements
    aux = 0;
    for i = 1:p
        aux = aux + a(i)*x(j+p-i);      % compute each signal's element by definition
    end
    x(j+p) = aux + noise(j+p);        % write it to the output vector (with noise)
end

x = x(p+1:end);   %cut the boundary conditions (which sort-of are not AR process)