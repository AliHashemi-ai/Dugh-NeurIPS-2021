function H = logdet(C)
%   LOGDET   Short description
%       [H] = LOGDET(C)
%
%   logdet is a computationally efficient operator that can deal with
%   sparse matrices


% assume diagonal form
%-----------------------------------------------------------------------
TOL   = 1e-8;                                        % c.f. n*max(s)*eps
n     = length(C);
s     = diag(C);
i     = find(s > TOL & s < 1/TOL);
C     = C(i,i);
H     = sum(log(diag(C)));

% invoke det if non-diagonal
%-----------------------------------------------------------------------
w     = warning;
warning off
[i j] = find(C);
if any(i ~= j)
      n = length(C);
      a = exp(H/n);
      H = H + log(det(C/a));
end
warning(w)

% invoke svd is rank deficient
%-----------------------------------------------------------------------
if imag(H) | isinf(H)
    s  = svd(full(C));
    H  = sum(log(s(s > TOL & s < 1/TOL)));
end

