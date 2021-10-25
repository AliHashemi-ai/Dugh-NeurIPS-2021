function [cov_mat,indice] = cov_mat_gen(dim,update_mode,variance)

%% 
% Five different modes for generating the spatial or temporal covariance matrices

% Input 
% 	dim : dimension of the squaure matrix (dim x dim)
% 	update mode: generate a specific model for the covariance matrix 
% 		update_mode = 'identity'; 
% 		update_mode = 'sparse'; 
% 		update_mode = 'full';
% 		update_mode = 'diagonal';
% 		update_mode = 'toeplitz';

% Output 
% 	cov_mat : covaricne matrix wiht identified dimension and structure 

%% 

switch update_mode
	case 'full'
		tmp = sqrt(variance) * randn(dim,2*dim);
		tmp = tmp*tmp';
		cov_mat = tmp/trace(tmp);
	case 'identity'	
		cov_mat = eye(dim);
	case 'diagonal'
		cov_mat = variance * diag(rand(dim,1));
	case 'sparse'
		N_active = 5;
		% % select active sources at random locations
		ind = randperm(dim);
		indice = ind(1:N_active);
		mat_cov_coef = zeros(dim,1); 
		mat_cov_coef(indice) = variance;
		cov_mat = spdiags(mat_cov_coef,0,dim,dim);  
	case 'toeplitz'
		% Generate the temporal correlation matrix with Toeplitz structure
		beta = 0.9;  % AR-coefficient (temporal corrolation parameter)
		c = beta.^abs([1:dim]'-1);
		r = beta.^abs(1-[1:dim]);
		cov_mat = toeplitz(c,r);
end

end