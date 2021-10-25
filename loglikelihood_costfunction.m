function [f] = loglikelihood_costfunction(X_total,Sigma)

% Calculating the value of the loglikelihood loss function for Type-II
% Bayesian learning

    [K,G] = size(X_total);
    f = log(det(Sigma))+(1/G)*sum(diag(X_total'*inv(Sigma)*X_total));
end