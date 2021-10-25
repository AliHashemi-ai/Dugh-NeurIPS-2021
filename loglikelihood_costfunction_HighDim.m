function [f] = loglikelihood_costfunction_HighDim(X_total, Gamma, B , Sigma)

% Calculating the value of the loglikelihood loss function for Type-II
% Bayesian learning for spatiotemporal setting. 

% Efficient implementation of "loglikelihood_costfunction" where the
% the Kronocker structure of covariance matrix Sigma_0 is simplified by 
% performing the matrix multiplication over its
% components B and Gamma. 

    [K, G] = size(X_total);
    N = size(Gamma,1);
    T = size(B,1);
    f = (T * log(det(Gamma))) + (N * log(det(B))) + (1/G)*sum(diag(X_total'*inv(Sigma)*X_total));
    % (Eq. 23 in the geodesic paper)
end