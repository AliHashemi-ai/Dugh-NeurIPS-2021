% Description:    This function estimates the covariance matrix of the data by imposing
%                 different kind of structures for spatial or temporal correlation matrices. 
	 
% Input

% Output:
% 	Sigma_structural: Final estimated covariance

function [Sigma_structural,p] = champ_colored_cov_est_realdata(C,M,Q,p,update_mode)

switch update_mode
    case 'Geodesic'
        % update_mode = 'Geodesic'; 
        S = inv(sqrtm(C))*sqrtm((sqrtm(C))*M*(sqrtm(C)))*inv(sqrtm(C));
%       S = sqrtm(C)*sqrtm(inv(sqrtm(C))*M*inv(sqrtm(C)))*sqrtm(C);
        S = S/trace(S);
        % (Eq. 18 in the ICLR paper)
        
%         % Efficient Implementation 
%         eps_default = 1e-8; 
%         [b_vec,b_val] = eig(C);
%         root_C_coeff = sqrt(max(real(diag(b_val)),0));
% 
%         inv_root_C_coeff = zeros(size(C,1),1);
%         inv_root_C_index = find(root_C_coeff >= eps_default);
%         inv_root_C_coeff(inv_root_C_index) = 1./root_C_coeff(inv_root_C_index);
% 
%         root_C = b_vec * diag(root_C_coeff) * b_vec';
%         inv_root_C = b_vec*diag(inv_root_C_coeff)*b_vec';
% 
%         [a_vec,a_val] = eig(root_C * M * root_C);
%         A_coeff = sqrt(max(real(diag(a_val)),0));
%         A = a_vec * diag(A_coeff) * a_vec';
%         S = inv_root_C * A * inv_root_C;
%         S = S/trace(S);            

    case 'Geodesic-notrace'
        % update_mode = 'Geodesic'; 
%       S = inv(sqrtm(C))*sqrtm((sqrtm(C))*M*(sqrtm(C)))*inv(sqrtm(C));
        % (Eq. 18 in the ICLR paper)
		
%       S = sqrtm(C)*sqrtm(inv(sqrtm(C))*M*inv(sqrtm(C)))*sqrtm(C);
%       S = S/trace(S);

        % Efficient Implementation 
        eps_default = 1e-8; 
        [b_vec,b_val] = eig(C);
        root_C_coeff = sqrt(max(real(diag(b_val)),0));

        inv_root_C_coeff = zeros(size(C,1),1);
        inv_root_C_index = find(root_C_coeff >= eps_default);
        inv_root_C_coeff(inv_root_C_index) = 1./root_C_coeff(inv_root_C_index);

        root_C = b_vec * diag(root_C_coeff) * b_vec';
        inv_root_C = b_vec*diag(inv_root_C_coeff)*b_vec';

        [a_vec,a_val] = eig(root_C * M * root_C);
        A_coeff = sqrt(max(real(diag(a_val)),0));
        A = a_vec * diag(A_coeff) * a_vec';
        S = inv_root_C * A * inv_root_C;
% 		
%       % S = S/trace(S);   uncomment for trace normalization constraint.      


%% ---------------------------------------------------------------------
%         % Efficient Implementation compatible with the paper
%         % Please change the C_M and C_N in the main code. Use their
%         % inverse, e.g., C_source = inv(L' * SigmaY_inv * L) and 
%         % C_noise = SigmaY_estimated. 
        
%         eps_default = 1e-8; 
%         [b_vec,b_val] = eig(C);
%         root_C_coeff = sqrt(max(real(diag(b_val)),0));
% 
%         inv_root_C_coeff = zeros(size(C,1),1);
%         inv_root_C_index = find(root_C_coeff >= eps_default);
%         inv_root_C_coeff(inv_root_C_index) = 1./root_C_coeff(inv_root_C_index);
% 
%         root_C = b_vec * diag(root_C_coeff) * b_vec';
%         inv_root_C = b_vec*diag(inv_root_C_coeff)*b_vec';
% 
% %         [a_vec,a_val] = eig(root_C * M * root_C);
%         [a_vec,a_val] = eig(inv_root_C * M * inv_root_C);
%         A_coeff = sqrt(max(real(diag(a_val)),0));
%         A = a_vec * diag(A_coeff) * a_vec';
%         S = root_C * A * root_C;
%% ---------------------------------------------------------------------        

 
    case 'Diagonal'
        % solving inner problem using diagonal constraint
        % (Eq. 20 in the ICLR paper)
        h = diag(C); 
        g = diag(M);    
        p = sqrt(g./(h));
        S = diag(p);

    case 'Toeplitz'
        % solving inner problem using circulant embedding 
        % (Eq. 53-55 in the geodesic paper)				
        h = diag(Q'*C*Q); % h is assigned based on the previous solution.     
        D = diag(Q'*C*M*C*Q);
        g = p.*D.*p;    
        ptmp = sqrt(g./(h));
%         S_tmp = Q*diag(ptmp)*Q';
%         ptmp = ptmp/trace(S_tmp);
        p = ptmp; 
        S = Q*diag(p)*Q';
end 
%%
Sigma_tmp = S ;
Sigma_structural = Sigma_tmp; 

end