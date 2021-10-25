function U = randorth(M)
% random orthogonal MxM matrix 
%
% Stefan Haufe, 2017

X = randn(M);
[Q,R] = qr(X);
R = diag(diag(R)./abs(diag(R)));
U = Q*R;