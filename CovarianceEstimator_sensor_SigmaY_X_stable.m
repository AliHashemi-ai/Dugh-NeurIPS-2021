% Description:    This function estimates the covariance matrix of the data by imposing
%                 different kind of structures for spatial or temporal correlation matrices. 
%				 

% Input
%     M:           size of SigmaY (Spatioal Correlation Matrix in sensor space)
%     T:           size of B (Temporal Correlation Matrix)
%     G:           Number of Block samples
%     Y:           NT-by-G data matrix
%    Sigma_0:      NT-by-NT matrx. True covariance. 

% Output:
% 	Sigma_structural: Final estimated covariance

% Reference: 
% Overleaf draft: NeurIPS paper% 

function [X_total,X_trial,Sigma_structural,Gamma,B,error_Gamma,error_B,f_loglikelihood]=CovarianceEstimator_sensor_SigmaY_X_stable(Y_total,L,SigmaY_0,Lambda_0,spatial_update_mode,temporal_update_mode)

N = size(L,2);
M = size(L,1);
G = size(Y_total,2);
T = size(Y_total,1) / M; 
% generate samples Xg

% Source Space 
% for g = 1:G
    % X(:,:,g)=reshape(X_total(:,g),T,N);
% end

block_length = T; 
% G = size(Y_total,2) / block_length; 
% % Number of Blocks
% %% Sensor Space 
% for g = 1:G
%     Y(:,:,g)=Y_total(:,((g-1)*block_length)+1:g*block_length)';
% end

%% Sensor Space
% G is the number of blocks. (it can be considered as different trials or different time windows making 
% the entire time frame.   
for g = 1:G
    Y(:,:,g)=reshape(Y_total(:,g),T,M);
end


%% ======================================
%%%% Initialize Gamma and B %%%%% =======
% =======================================
%% Pseudo-Inverse Initialization
% 
% L_square=sum(L.^2,1);              
% Inv_L_square=zeros(1,N);
% Index_pos=find(L_square>0);
% Inv_L_square(Index_pos)=1./L_square(Index_pos);
% w_filter=spdiags(Inv_L_square',0,N,N)*L';
% gamma_init=mean(mean((w_filter*Y_total).^2));
% gammas=gamma_init*ones(N,1);
% Gamma_init = spdiags(gammas,0,N,N);
% Gamma_init = Gamma_init/trace(Gamma_init);

%% Random Initialization
% tmp = rand(N,2*N);
% Gamma_init = tmp*tmp';
% Gamma_init = Gamma_init/trace(Gamma_init);
% X_init = randn(N,size(Y_total,2));
% X_total = X_init;

%% Initializaation with least square (\ell_2) psudo-inverse solution. 
% Spatial matched-fileter
L_sqaure = sum(L.^2,1);
inv_L_sqaure = zeros(1,N); 
L_nonzero_index = find(L_sqaure > 0); 
inv_L_sqaure(L_nonzero_index) = 1./L_sqaure(L_nonzero_index);
w_filter = spdiags(inv_L_sqaure',0,N,N)*L';
Y_ave = mean(Y,3)';
gammas = mean(mean(w_filter * Y_ave) .^2);
Gamma_init = double(gammas)*speye(N,N);
% Gamma_init = speye(N,N);
Gamma_k = Gamma_init;     
gamma_k = gammas; 
% Cov_noise = stdnoise^2 * eye(M); % diagonal and scalar
% % Cov_noise = stdnoise * ones(M,1);
% % Cov_noise = spdiags(stdnoise',0,M,M);

% SigmaY = Cov_noise + L*Gamma_k*L'; % Covariance of Y
% SigmaY_init = SigmaY; 
% % Invert CM keeping symmetry
% [U,S,V] = svd(SigmaY);
% SigmaY_inv = U*diag(1./(diag(S)+eps))*U';		 
% SigmaY_invL = SigmaY_inv * L;			
% Linv = Gamma_k * L' * SigmaY_inv;
% X_total = Linv * Y_total;
% X_old = X_total;     

% tmp = rand(N,100);
% Gamma_init = tmp*tmp';
% SigmaY_init = (stdnoise^2 * eye(M)) + (L * Gamma_init *L');
% SigmaY_init = SigmaY_init/trace(SigmaY_init);

% Sigma_Y_initialization 
update_mode = 1; 
if update_mode == 1
    SigmaY_init = Lambda_0 + (L * Gamma_init *L');
%     SigmaY_init = (stdnoise^2 * eye(M)) + (L * Gamma_init *L');
else
    tmp = rand(M,2*M);
    SigmaY_init = tmp*tmp';
    SigmaY_init = SigmaY_init/trace(SigmaY_init);
end
%% --- Temporal matrix initialization --- %% 
%% Toeplitz Initialization
% beta = 0.8;
% c = beta.^abs([1:T]'-1);
% r = beta.^abs(1-[1:T]);
% B_init = toeplitz(c,r);
% B_init = B_init/trace(B_init);

%% Random Initialization 
tmp =  rand(T,T*2);
B_init = tmp*tmp';
B_init = B_init/trace(B_init);

%% Identity Initialization 
% B_init = eye(T);

%% ======================================================== %% 
% initialization for the temporal setting.
% Use Toeplits structure for diagonalization by Fourier transform.
L_embedding = 2*T-1;
Q_time = [eye(T) zeros(T,L_embedding -T)]*dftmtx(L_embedding )/sqrt(L_embedding);
% B_inv = inv(B_k);
% h = diag(Q_time'*B_inv*Q_time);           
p_int = rand(L_embedding -1,1);
p_int = (p_int+flipud(p_int))/2;
p_int = [rand(1); p_int];
p_time = p_int;

% initialization for the spatial setting.
stdnoise = mean(diag(Lambda_0));
if stdnoise ~= 0 
    Q_space = [L eye(M)];
    %% initialize the signal variances %% 
%   p_space = rand(N+M,1); % random 
    p_space = ones(N+M,1); % with ones 
%   p_space = vec_init *ones(N+M,1);    % with Psudo-Invers 

    %% initialize the noise variances %% 
	p_space(N+1:end,1) = diag(Lambda_0);
%    p_space(N+1:end,1) = (stdnoise^2) * ones(M,1);
    % p_space(N+1:end,1) = stdnoise' .* ones(M,1);
else
    Q_space = L;
    p_space = rand(N,1);
end
H_k = diag(p_space); 

%% ======================================================== %% 

Sigma_tmp = kron(SigmaY_init,B_init);
Sigma_tmp = Sigma_tmp / norm(Sigma_tmp, 'fro'); 
% for g = 1:G
%     Y(:,:,g)=reshape(Y_total(:,g),M,T);
% end
% Y_loss = reshape(Y,M,[]);
f_loglikelihood =  loglikelihood_costfunction_HighDim(Y_total,SigmaY_init, B_init, Sigma_tmp); 

SigmaY_k = SigmaY_init;
B_k = B_init;
% P_k = diag(p_time);
Converge_X = false;
EPSILON = 1e-5; 
Max_num_iterations = 700; 
itr_count = 0;
PRINT = 1; 
X_old = 0;
gamma_old = 0;

while(~Converge_X)
	% Setting the threshold based on NMSE of estimating the temporal and
	% spatial correlation matrices
	tol = 1e-2;
	Converge = 0;

	while(~Converge)
		% fix SigmaY update B
		M_time=0;
		%% Sensor Space
		for g = 1:G
			M_time =M_time+Y(:,:,g)*inv(SigmaY_k)*Y(:,:,g)';
		end
% 		M_time = M_time/G*T;  
		M_time = M_time/G*M; % (Eq. 37 in the geodesic paper)
		% temporal_update_mode = 'Geodesic';
		
		switch temporal_update_mode
		
			case 'Geodesic'
				% Update B by finding the geometric mean between the whitened version of sample convergence matrix and previous solution of B
				B = sqrtm(B_k)*sqrtm(inv(sqrtm(B_k))* M_time *inv(sqrtm(B_k)))*sqrtm(B_k);
                B = B / trace(B);
				% (Eq. 9 in the NeurIPS paper)				
                
% %                 D = diag(Q_time'*inv(B_k)*M_time*inv(B_k)*Q_time);					
% %                 M_time_eff = P_k.*D.*P_k;    
%                 
%               D = P_k*Q_time'*inv(B_k)*M_time*inv(B_k)*Q_time*P_k';					
% 				P = sqrtm(P_k)*sqrtm( inv(sqrtm(P_k))* D *inv(sqrtm(P_k)) )*sqrtm(P_k);
% 				B = Q_time*P*Q_time';

			
			case 'Toeplitz'
				% solving inner problem using circulant embedding 
				% (Eq. 19-21 in the NeurIPS paper)
				
				% L_embedding = 2*T-1;
				% Q = [eye(T) zeros(T,L_embedding -T)]*dftmtx(L_embedding )/sqrt(L_embedding);
				% B_inv = inv(B_k);
				% h = diag(Q'*B_inv*Q);           
				% p_int = rand(L_embedding -1,1);
				% p_int = (p_int+flipud(p_int))/2;
				% p_int = [rand(1); p_int];
				% p = p_int;
				
                Converge_p_time = 0;
                itr_p_time = 0;
                tol_B = 1e-3;
                
                B_inv = inv(B_k);
                z = diag(Q_time'*B_inv*Q_time); % z is assigned based on the previous solution.     
				
				% B_inv = inv(B_k);
%                 z = diag(P_k); % h is assigned based on the previous solution.
				
				while(~Converge_p_time)
					itr_p_time = itr_p_time+1;  
% % 					B_inv = inv(B_k);
% % 					z = diag(Q_time'*B_inv*Q_time);      
				
					B = Q_time*diag(p_time)*Q_time';
					D = diag(Q_time'*inv(B)*M_time*inv(B)*Q_time);
%                     D = diag(Q_time'*inv(B)'*M_time*inv(B)*Q_time);
					
					g = p_time.*D.*p_time;    
					ptmp = sqrt(g./(z));
					
					
% 					D = diag(Q_time'*inv(B_k)*M_time*inv(B_k)*Q_time);					
% 					g = p_time.*D.*p_time;    
% 					ptmp = sqrt(g.*(z));
					
					B_tmp = Q_time*diag(ptmp)*Q_time';
					ptmp = ptmp/trace(B_tmp);
					Converge_p_time = (norm(ptmp-p_time)<tol_B);
                    
					if mod(itr_p_time,10) == 0 
						fprintf('p_time_error: %d\n',(norm(ptmp-p_time)));
					end
					
					p_time = ptmp; 
				end
				B = Q_time*diag(p_time)*Q_time';
%               P = diag(p_time);
		end 
	%   Converge = norm(B-B_k,'fro')<tol;
	
		fprintf('B_error: %d \n \n',norm(B-B_k,'fro'));
        error_B(itr_count+1) = norm(B-B_k,'fro'); 

		Converge_B = norm(B-B_k,'fro')<tol;
		
	%   B_k = B/trace(B); % Normalize it based on its trace. 
		B_k = B;
%       P_k = P;
		% Sigma_tmp = kron(SigmaY_k,B_k);

	%%  fix B update Gamma and SigmaY
		
		M_space=0;
		for g = 1:G
			M_space =M_space+Y(:,:,g)'*inv(B_k)*Y(:,:,g);
		end
% 		M_space = M_space/G*M;
        M_space = M_space/G*T;
	%   M_space = M_space/G;
	%   M_space = M_space/G*N; %% In case of updating rules in the source space 	
		switch spatial_update_mode
			case 'Geodesic'
				% spatial_update_mode = 'Geodesic'; 
% 				SigmaY = sqrtm(SigmaY_k)*sqrtm(inv(sqrtm(SigmaY_k))*M_space*inv(sqrtm(SigmaY_k)))*sqrtm(SigmaY_k);

                eps_default = 1e-8;
                SigmaY = Q_space * H_k * Q_space';
                [U,S,V] = svd(SigmaY);
                SigmaY_inv = U*diag(1./(diag(S)+eps_default))*U';
                M_SN = H_k' * Q_space' * SigmaY_inv' * M_space * SigmaY_inv * Q_space * H_k;	
%               M_SN = H_k * Q_space' * SigmaY_inv * M_space * SigmaY_inv * Q_space * H_k;	

%                 SigmaY = Q_space * H_k * Q_space';
%                 [U,S,V] = svd(SigmaY);
%                 SigmaY_inv = U*diag(1./(diag(S)+eps))*U';
%                 L_inv = H_k * Q_space' * SigmaY_inv;
%                 for g = 1:G
%                     mu(:,((g-1)*block_length)+1:g*block_length)= L_inv *  Y(:,:,g)' * sqrtm(inv(B_k));
%                     X_trial(:,:,g) = L_inv *  Y(:,:,g)' * sqrtm(inv(B_k));
%                 end
% %               X_total = L_inv * Y_total;    
%                 X_total = mu;  
% 
%                 M_SN_test1 = mean(X_total.^2,2);
%                 M_SN_test2 = sum(L_inv * M_space .* L_inv, 2);
  
				% Efficient & more stable
		        eps_default = 1e-8; 
		        [b_vec,b_val] = eig(H_k);
		        root_C_coeff = sqrt(max(real(diag(b_val)),0));
		
		        inv_root_C_coeff = zeros(size(H_k,1),1);
		        inv_root_C_index = find(root_C_coeff >= eps_default);
		        inv_root_C_coeff(inv_root_C_index) = 1./root_C_coeff(inv_root_C_index);
		
		        root_C = b_vec * diag(root_C_coeff) * b_vec';
		        inv_root_C = b_vec*diag(inv_root_C_coeff)*b_vec';
		
		        [a_vec,a_val] = eig(inv_root_C * M_SN * inv_root_C);
		        A_coeff = sqrt(max(real(diag(a_val)),0));
		        A = a_vec * diag(A_coeff) * a_vec';
		        S = root_C * A * root_C;
                H = S; 
                
%               %%% sqrtm function
				H = sqrtm(H_k)* sqrtm(inv(sqrtm(H_k))*M_SN*inv(sqrtm(H_k))) * sqrtm(H_k);
             

%               H = H / trace(H);
                H_k = H;
                % (Eq. 12 in the NeurIPS paper)

%               Gamma = real(H(1:N,1:N));
                Gamma = H(1:N,1:N);
                Converge_Gamma = norm(Gamma-Gamma_k,'fro')<tol;
                error_Gamma(itr_count + 1) = norm(Gamma-Gamma_k,'fro'); 
                fprintf('Gamma_error: %d \n \n ',norm(Gamma-Gamma_k,'fro'));
                Gamma_k = Gamma;
                gammas = (diag(Gamma_k));
                
				
			case 'sum_of_rank_one'
				% solving inner problem using sum of rank one constraint
				% This is an updating rule based on the augmented approach that
				% we developed based on the MM framework. Please check the
				% updating rules for gamma in MM paper. (Appendix D)
				
				% itr_p = 0;	
				% if stdnoise ~= 0 
					% Q = [L eye(M)];
					% p = rand(N+M,1);
					% p(N+1:end,1) = (stdnoise^2) * ones(M,1);
				% else
					% Q = L;
					% p = rand(N,1);
				% end

				SigmaY_inv = inv(SigmaY_k);       
	%           [U,S,V] = svd(SigmaY_k);
	%			SigmaY_inv = U*diag(1./(diag(S)+eps))*U';

	% 			f_previous = loglikelihood_costfunction_HighDim(Y_total,SigmaY_k, B_k, kron(SigmaY,B_k));

    
				check_idea = false;
%                 check_idea = true;

                homoscedastic_update = false;
                
%               H_k = diag(p_space);
				if check_idea
					%% Check the idea of full-structure covariance                 
					z = p_space;
				else 	
					z = diag(Q_space'*SigmaY_inv*Q_space);
				end

					
				itr_p = 0;
				Converge_p_space = 0;
				tol_p_SigmaY = 1e-2;
				
				while(~Converge_p_space)
					itr_p = itr_p+1;
					SigmaY = Q_space * diag(p_space) *Q_space';
%                   SigmaY = Q_space*H_k*Q_space';
% % 				f_previous = loglikelihood_costfunction_HighDim(Y_total,SigmaY_k, B_k, kron(SigmaY,B_k));
					D = diag(Q_space'* inv(SigmaY) * M_space * inv(SigmaY) * Q_space);
					g = p_space.*D.*p_space;
					
					if check_idea
						%% Check the idea of full-structural covariance                 
						ptmp = sqrt(g.*z);
                    else 					
						if homoscedastic_update 
							ptmp(1:N,1) = (sqrt(g(1:N)./z(1:N)));
							ptmp(N+1:N+M,1) = ones(M,1)*sqrt(sum(g(N+1:N+M))/sum(z(N+1:N+M)));
						else
							ptmp = sqrt(g./z);
						end  
					end
					
%                   ptmp = sqrt(g./h);
% 					SigmaY_tmp = Q_space*diag(ptmp)*Q_space';
% 					ptmp = ptmp/trace(SigmaY_tmp);
					Converge_p_space = (norm(ptmp-p_space)<tol_p_SigmaY);		
	%                 fprintf('p_error: %d\n',(norm(ptmp-p)));
					if mod(itr_p,20) == 0 
						fprintf('p_error: %d\n',(norm(ptmp-p_space)));
					end
	% 				f_current = loglikelihood_costfunction(Y_total,kron(SigmaY_tmp,B_k));
	% 				Converge = abs(f_current-f_previous)/max(1,abs(f_previous))<tol;
					p_space = ptmp;
				end
				SigmaY = Q_space*diag(p_space)*Q_space';
%               SigmaY = Q_space*p_space*Q_space';
%               SigmaY = Q_space*H_k*Q_space';
                gammas = (real(p_space(1:N)));
                
                Converge_Gamma = max(abs(gammas - gamma_k)) <tol;
                error_Gamma(itr_count + 1) = max(abs(gammas - gamma_k)); 
                fprintf('gamma_error:%d\n\n ',max(abs(gammas - gamma_k)));
        		gamma_k = gammas;
		end 
	%================== Checking convergence for SigmaY ======================%% 
	%   Converge = norm(SigmaY-SigmaY_k,'fro')<tol;

	%================== Calculating the error in updating SigmaY ===============
		Converge_SigmaY = norm(SigmaY-SigmaY_k,'fro')<tol;
        error_SigmaY(itr_count + 1) = norm(SigmaY-SigmaY_k,'fro'); 
		fprintf('SigmaY_error:%d\n\n ',norm(SigmaY-SigmaY_k,'fro'));

    %================== Calculating the error in updating Gamma ===============


        
	%====================== Update the SigmaY and Gamma ==================================
		%%SigmaY_k = SigmaY/trace(SigmaY);
		SigmaY_k = SigmaY;


	%% Convergence Analysis 

	% Convergence based on the estimation error in spatioal or temporal correlation
        Converge =  (Converge_Gamma) ||  (Converge_B);
% 		Converge =  (Converge_SigmaY) ||  (Converge_B);
%       Converge =  (Converge_SigmaY) &&  (Converge_B); % Converge slowly  
		
	% =========================================================================
	% Convergence based on the negative loglikelihood cost function value  
		% Sigma_tmp = kron(SigmaY_k, B_k);
		% f_loglikelihood  = [f_loglikelihood, loglikelihood_costfunction_HighDim(Y_total,SigmaY, B, Sigma_tmp)];
		
        Sigma_tmp = kron(SigmaY_k, B_k);
		f_loglikelihood  = [f_loglikelihood, loglikelihood_costfunction_HighDim(Y_total,SigmaY, B, Sigma_tmp)];

        
        % Converge = abs(f_loglikelihood(end)-f_loglikelihood(end-1))/max(1,abs(f_loglikelihood(end-1)))<tol;
	%   fprintf('convergence_error: %d \n \n ',abs(abs(f_loglikelihood(end))-abs(f_loglikelihood(end-1)))/max(1,abs(f_loglikelihood(end-1))));
	% %   Converge = abs(abs(f_loglikelihood(end))-abs(f_loglikelihood(end-1)))/max(1,abs(f_loglikelihood(end-1)))<tol;

	% =========================================================================
	% Convergence based on the value of the ground-truth
	% %   f_previous = loglikelihood_costfunction(Y_total,Sigma_0);
	% %   f_current = loglikelihood_costfunction(Y_total,Sigma_tmp);
	% %   Converge = abs(f_current-f_previous)/max(1,abs(f_previous))<tol;

	end
	
% 	X_update_mode = 'Approx_stable';
    X_update_mode = 'Exact';
	switch X_update_mode 
		case 'Approx_classic'
			Gamma_tilda = diag(p);
			Gamma = Gamma_tilda(1:N,1:N);
			SigmaY_estimated = (L * Gamma * L') + Lambda_0;
			X_estimated = Gamma * L' * inv(SigmaY_estimated) * Y_total;   
			X_total = X_estimated; 
			
		case 'Approx_stable'
			% Invert CM keeping symmetry
            block_length = T;
			gammas_diag = spdiags(p_space(1:N),0,N,N); 
            % Gamma = spdiags(p(1:N),0,N,N); 
            % SigmaY_estimated = (L * Gamma * L') + ((stdnoise^2 * eye(M)));
            % SigmaY = SigmaY_estimated; 
            % SigmaY_k = SigmaY_estimated; 
			[U,S,V] = svd(SigmaY);
			SigmaY_inv = U*diag(1./(diag(S)+eps))*U';
% 			SigmaY_invL = SigmaY_inv * L;

            switch spatial_update_mode
                case 'Geodesic'
                      L_inv = Gamma * L' * SigmaY_inv;
%                     L_inv = H * Q_space' * SigmaY_inv;
                case 'sum_of_rank_one' 
                    L_inv = gammas_diag * L' * SigmaY_inv;
            end
            
            for g = 1:G
%                 mu(:,((g-1)*block_length)+1:g*block_length)= L_inv * Y(:,:,g)' ;
%                 X_trial(:,:,g) = L_inv * Y(:,:,g)';
                
                mu(:,((g-1)*block_length)+1:g*block_length)= L_inv *  Y(:,:,g)' * sqrtm(inv(B_k));
                X_trial(:,:,g) = L_inv *  Y(:,:,g)' * sqrtm(inv(B_k));
                
            end
% 			X_total = L_inv * Y_total;    
			X_total = mu;  
%           X_total = X_trial;  
%           X_total = mean(X_trial,3);  
			
			% X_estimated_SigmaY = pinv(L) * (SigmaY_estimated - (stdnoise^2 * eye(M)) ) * inv(SigmaY_estimated) * Y_total;
		
		case 'Exact'
		% 	Check Section 3.3, "Efficient Solution for the Posterior Mean" of the Geodesic Paper in the overleaf: "Elsevier-NeuroImage-Journal II-Main"
		%   lambda = (mean(p(N+1:end)));
			lambda = p_space(N+1:N+M,1);
            
            Gamma = spdiags(p_space(1:N),0,N,N); 
			P = zeros(T,M);
			mu_eff = zeros(T,N);
			[Vx,lambda_x] = eig(L*Gamma*L');
			[Vt,lambda_t] = eig(B_k*eye(T)); 
			lambda_x = spdiags(diag(lambda_x),0,length(lambda_x),length(lambda_x));
			lambda_t = spdiags(diag(lambda_t),0,length(lambda_t),length(lambda_t));

			for  Eff_index_i = 1:M
					for Eff_index_j =1:T
					%   P(Eff_index_i,Eff_index_j) = 1/(lambda_x(Eff_index_i,Eff_index_i) * lambda_t(Eff_index_j,Eff_index_j) + 1); 
						temp = ( (lambda_x(Eff_index_i,Eff_index_i) * lambda_t(Eff_index_j,Eff_index_j)) + lambda(Eff_index_i)); 
						P(Eff_index_j,Eff_index_i) = 1/ temp;              
					end
			end

% 			if block_length == size(Y_total,2)
% 				mu_eff =  Vt * lambda_t* (P .* (Vt'* Y_total' * Vx)) * Vx' * L * Gamma';
% 			else 
			%% Sensor Space 
				for g = 1:G
					mu_eff(((g-1)*block_length)+1:g*block_length,:)= Vt * lambda_t* (P .* (Vt'* Y(:,:,g) * Vx)) * Vx' * L * Gamma';
                    X_trial(:,:,g) = Vt * lambda_t* (P .* (Vt'* Y(:,:,g) * Vx)) * Vx' * L * Gamma';
				end
% 			end
			% mu = reshape(mu_eff,N,T);
			X_estimated = mu_eff';  
			X_total = X_estimated; 
	end
	
	
	% ================= Check stopping conditions, etc. ============== 
	
	itr_count = itr_count + 1;
 
	usedNum = nnz(p_space(1:N)); 
% 	gammas = p_space(1:N); 
    if (PRINT) disp(['iters: ',num2str(itr_count),'   num coeffs: ',num2str(usedNum), ...
            '   gamma change: ',num2str(max(abs(gammas - gamma_old)))]); 
    end

    if (size(X_total) == size(X_old))
		fprintf('dX_Frob_norm: %d \n \n',norm(X_old-X_total,'fro'));
		Converge_X = norm(X_old-X_total,'fro') < EPSILON;
		dX_new = max(max(abs(X_old - X_total)));
		fprintf('dX_Abs_norm: %d \n \n',dX_new);
        if (PRINT==1) 
%           disp(['   mu change: ',num2str(max(max(abs(mu_old - mu))))]);
            error_X = norm(X_total-X_old,'fro')^2/norm(X_old,'fro')^2;
            disp(['   X change: ',num2str(max(max(abs(error_X))))]);
        end
        if dX_new < EPSILON || (itr_count > Max_num_iterations)
            break;  
			Converge_X = true; 
        end
    end
	
	X_old = X_total; 
	gamma_old = gammas;
end


%%
Gamma = diag(real(p_space(1:N)));
B = B_k; 
Sigma_tmp = kron(SigmaY_k, B_k);
Sigma_structural = Sigma_tmp; 

    % figure
    % plot(1:length(f_loglikelihood),abs(f_loglikelihood),'b-');
    % xlabel('Iteration');
    % ylabel('Objective Value');
    % title('The convergence behaviour tracking');

end