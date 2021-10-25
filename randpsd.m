function [sqC, C, eigvals, Q, condi] = randpsd(M, r, alpha)
% random MxM positive semidefinite matrix 
%
% r: number of nonzero eigenvalues, default: M
% alpha: range of eigenvalue scales, determines condition if r=M
%
% sqC: decomposition such that sqC*sqC' = C
% C: random PSD matrix
% eigvals: nonzero eigenvalues of C
% Q: correspondong columns of random orthogonal matrix
% condi: condition of C
%
% Stefan Haufe, 2017

if nargin < 2
  r = M;
end

if nargin < 3
  alpha = 2;
end

Q = randorth(M);
Q = Q(:, 1:r);
eigvals = exp(alpha*randn(r, 1));
eigvals = eigvals / max(eigvals);

sqC = diag(sqrt(eigvals))*Q';

if nargout > 1
  C = sqC'*sqC;
end
% sqC = sqrtm(C);

if r == M
  condi = (max(eigvals)/min(eigvals));
else
  condi = inf;
end

