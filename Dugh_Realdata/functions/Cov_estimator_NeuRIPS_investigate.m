function [X,gamma_est,Energy,Lambda,W] = Cov_estimator_NeuRIPS_investigate(Y_total,L, varargin)
% Champagne algorithm for hetrosedasdic noise setting. 

%% Sensor Space
% G is the number of blocks. (it can be considered as different trials or different time windows making 
% the entire time frame.   

    

% Dimension of the Problem
D = 1; 
[M,ND] = size(L); 
N = ND /D;  

G = size(Y_total,2);
T = size(Y_total,1) / M; 

block_length = T; 
for g = 1:G
    Y_block(:,:,g)=reshape(Y_total(:,g),T,M);
    Y(:,((g-1)*block_length)+1:g*block_length) = Y_block(:,:,g)';
end


% Default Control Parameters  
EPSILON = 1e-8;       			% threshold for stopping iteration. 
Max_num_iterations = 50;        % maximum iterations
PRINT = 0;          			% don't show progress information
print_figures = 1; 				% Plot the amplitude of the sources and their corresponding estimation
								% along with the plots for oroginal and estimated noise covariances
itr_count = 1;
noise_update_mode = 'Diagonal'; 
source_update_mode = 'Diagonal'; 
% source_update_mode = 'Geodesic-notrace'; 
% source_update_mode = 'Geodesic'; 

% X_init = w_filter*Y; 
% E_init = Y - L * X_init; 
% Lambda_init = cov(E_init');

% Peloreta = mkfilt_eloreta_v2(L,0.05);    
% Weight_eloreta = Peloreta'*Y; 
% E_init_elor = Y - L * Weight_eloreta; 
% Lambda_init = cov(E_init_elor');
% gamma_init=mean(mean((Weight_eloreta*Y).^2));
% gammas = gamma_init*ones(N,1);
% Gamma_init = spdiags(gammas,0,N,N);

YYt = Y * Y'; 
C_y = 1/T * YYt; % Sample covarinace matrix  

tmp = rand(M,T);
Lambda_init = tmp*tmp';
% Lambda_init = eye(size(C_y))*max(eig(C_y))*1e-3;
% Lambda_init = Lambda_init/trace(Lambda_init);

	
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
 
for i=1:2:(length(varargin)-1)
    switch lower(varargin{i})
        case 'noise_cov'
            Lambda_init = varargin{i+1}; 
		case 'noise_update_mode'
			noise_update_mode = varargin{i+1};	
        case 'temporal_update_mode'
            temporal_update_mode = varargin{i+1};	
        case 'epsilon'   
            EPSILON = varargin{i+1}; 
        case 'print'    
            PRINT = varargin{i+1}; 
        case 'print_figures'    
            print_figures = varargin{i+1}; 
        case 'max_iters'
            Max_num_iterations = varargin{i+1};  
        otherwise
            error(['Unrecognized parameter: ''' varargin{i} '''']);
    end
end

%% ======================================================== %% 
% initialization for the spatial setting.
% Q_source = eye(N);
% p_source = rand(N,1);

Q_source = eye(ND);
p_source = rand(ND,1);

Q_noise = eye(M);
p_noise = rand(M,1);

if strcmp(noise_update_mode,'Toeplitz')
    L_embedding = 2*M-1;
    Q_noise = [eye(M) zeros(M,L_embedding-M)]*dftmtx(L_embedding)/sqrt(L_embedding);      
    p_noise = rand(L_embedding -1,1);
    p_noise = (p_noise+flipud(p_noise))/2;
    p_noise = [rand(1); p_noise];
end
 
%% ===========================================
%%%% Initialize Lambda and Gamma %%%%% =======
% ============================================

%% Random Initialization 
% tmp = rand(N,2*N);
% Gamma_init = tmp*tmp';
% Gamma_init = Gamma_init/trace(Gamma_init);

%% Pseudo-Inverse Initialization for Gamma
L_square=sum(L.^2,1);              
Inv_L_square=zeros(1,ND);
Index_pos=find(L_square>0);
Inv_L_square(Index_pos)=1./L_square(Index_pos);
w_filter=spdiags(Inv_L_square',0,ND,ND)*L';
% X_init = w_filter * Y; 
% Gamma_init = X_init * X_init'/ T; 
gamma_init=mean(mean((w_filter*Y).^2));
gammas=gamma_init*ones(ND,1);
% p_source = gammas; 
Gamma_init = spdiags(gammas,0,ND,ND);
% % Gamma_init = Gamma_init/trace(Gamma_init);

Lambda = Lambda_init; 
Gamma = Gamma_init; 


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

X_old = 0;
print_convergence = false;
Converge_X = false;   

keep_list = [1:ND]';
active_set = keep_list'; 
usedNum = length(keep_list);
B_k = B_init;

% *** Learning loop ***
while(~Converge_X)

% 	pruning_threshold = -1; % Active-set strategy
%     % *** Prune weights as their hyperparameters go to zero ***
%     if (min(p_source) < pruning_threshold)
%         index = find(p_source > pruning_threshold);
%         p_source = p_source(index);  % use all the elements larger than MIN_GAMMA to form new 'gamma'
%         L = L(:,index);    % corresponding columns in Phi
%         keep_list = keep_list(index);
%         usedNum = length(p_source);
%     end 
	
	Gamma = spdiags(p_source,0,usedNum,usedNum);
	SigmaY_estimated = (L * Gamma * L') + (Lambda);
    

%     eps_default = 1e-8; 
%     [b_vec,b_val] = eig(SigmaY_estimated);
%     root_C_coeff = sqrt(max(real(diag(b_val)),0));
% 
%     inv_root_C_coeff = zeros(size(SigmaY_estimated,1),1);
%     inv_root_C_index = find(root_C_coeff >= eps_default);
%     inv_root_C_coeff(inv_root_C_index) = 1./root_C_coeff(inv_root_C_index);
%     inv_root_C = b_vec*diag(inv_root_C_coeff)*b_vec';
%     SigmaY_inv = inv_root_C; 
    
	SigmaY_inv = inv(SigmaY_estimated);  
	X_estimated = Gamma * L' * SigmaY_inv * Y;   
	X_total = X_estimated;    
	C_source = L' * SigmaY_inv * L; 

%     
% 	M_source=0;
% 	M_source = M_source + X_total * X_total';
% 	M_source = M_source/T;	 % (Eq. 11 in the ICLR paper)
	
    Gamma_old = Gamma; 
% 	[Gamma,p_source] = champ_colored_cov_est_realdata(C_source,M_source,Q_source,p_source,source_update_mode); 
    
%     X_B = X_total * inv(B_k);
%     x2=mean(X_B .* X_total,2);
    x2=mean(X_total.^2,2);
    fc = L' * SigmaY_inv;  
    z=sum(fc.*L',2);
    vcs = 1;
    if vcs == 0
        x20=sum(reshape(x2,D,N),1);
        z0=sum(reshape(z,D,N),1);
        v0=zeros(size(z0));
        ff=find(z0>0);
        v0(ff)=sqrt(max(x20(ff)./z0(ff),0)); % CN 10/2012 added max,0
    %   v0(ff)=sqrt(x20(ff)./z0(ff));
        p_source=reshape(ones(D,1)*v0,ND,1);
    else
%       vvec=zeros(size(x2));
        ff=find(z>0);
        p_source(ff)=sqrt(max(x2(ff)./z(ff),0)); % CN 10/2012 power 4.88 update v
    %   vvec(ff)=sqrt(x2(ff)./z(ff));
    end    

    if print_convergence 
        fprintf('Gamma_error: %d \n \n ',norm(Gamma-Gamma_old,'fro'));
    end
    
    Gamma = spdiags(p_source,0,usedNum,usedNum);
    SigmaY_estimated = (L * Gamma * L') + (Lambda);
    [p,d]=eig(SigmaY_estimated);
    d=max(real(diag(d)),0);
    
%     [b_vec,b_val] = eig(SigmaY_estimated);
%     root_C_coeff = sqrt(max(real(diag(b_val)),0));
% 
%     inv_root_C_coeff = zeros(size(SigmaY_estimated,1),1);
%     inv_root_C_index = find(root_C_coeff >= eps_default);
%     inv_root_C_coeff(inv_root_C_index) = 1./root_C_coeff(inv_root_C_index);
%     inv_root_C = b_vec*diag(inv_root_C_coeff)*b_vec';
%     SigmaY_inv = inv_root_C; 
    
	SigmaY_inv = inv(SigmaY_estimated);  
	X_estimated = Gamma * L' * SigmaY_inv * Y;   
	X_total = X_estimated;    
	E_total = Y - L * X_total; 
	M_noise = 0;
	M_noise = M_noise + E_total * E_total';
	M_noise = M_noise/T;	 % (Eq. 18 in the ICLR paper)

    Lambda_old = Lambda;
	C_noise = SigmaY_inv; 
%   C_noise = SigmaY_estimated; 
	[Lambda,p_noise] = champ_colored_cov_est_realdata(C_noise,M_noise,Q_noise,p_noise,noise_update_mode);

    if print_convergence 
        fprintf('Lambda_error: %d \n \n ',norm(Lambda-Lambda_old,'fro'));
    end
     
    SigmaY_k = (L * Gamma * L') + (Lambda);
    SigmaY = SigmaY_k;
    % fix SigmaY update B
    M_time=0;
    %% Sensor Space
    for g = 1:G
        M_time =M_time+Y_block(:,:,g)*inv(SigmaY_k)*Y_block(:,:,g)';
    end
% 		M_time = M_time/G*T;  
    M_time = M_time/G*M; % (Eq. 37 in the geodesic paper)
    % temporal_update_mode = 'Geodesic';

    switch temporal_update_mode
        
        case 'Identity' 
            B = eye(T);
            
        case 'Geodesic'
            % Update B by finding the geometric mean between the whitened version of sample convergence matrix and previous solution of B
            B = sqrtm(B_k)*sqrtm(inv(sqrtm(B_k))* M_time *inv(sqrtm(B_k)))*sqrtm(B_k);
            B = B / trace(B);
            % (Eq. 9 in the NeurIPS paper)				

        case 'Toeplitz'
            % solving inner problem using circulant embedding 
            % (Eq. 19-21 in the NeurIPS paper)

            Converge_p_time = 0;
            itr_p_time = 0;
            tol_B = 1e-1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE THIS

%           B_inv = inv(B_k);
            [U,S,V] = svd(B_k);
            B_inv = U*diag(1./(diag(S)+eps))*U';
            z = diag(Q_time'*B_inv*Q_time); % z is assigned based on the previous solution.     

            while(~Converge_p_time)
                itr_p_time = itr_p_time+1;  
% % 			B_inv = inv(B_k);

%               [U,S,V] = svd(B_k);
%               B_inv = U*diag(1./(diag(S)+eps))*U';
            
% % 			z = diag(Q_time'*B_inv*Q_time);      

                B = Q_time*diag(p_time)*Q_time';
                [U,S,V] = svd(B);
                B_inv = U*diag(1./(diag(S)+eps))*U';
                Diag_M_time = diag(Q_time'*B_inv*M_time*B_inv*Q_time);
                
%                 Diag_M_time = diag(Q_time'*inv(B)*M_time*inv(B)*Q_time);

                g = p_time.*Diag_M_time.*p_time;    
                ptmp = sqrt(g./(z));

                B_tmp = Q_time*diag(ptmp)*Q_time';
                ptmp = ptmp/trace(B_tmp);
                Converge_p_time = (norm(ptmp-p_time)<tol_B);
                
                if PRINT == 1
                    if mod(itr_p_time,10) == 0 
                        fprintf('p_time_error: %d\n',(norm(ptmp-p_time)));
                    end
                end

                p_time = ptmp; 
            end
            B = Q_time*diag(p_time)*Q_time';
%               P = diag(p_time);
    end 
    if PRINT == 1
        fprintf('B_error: %d \n \n',norm(B-B_k,'fro'));
    end 
%   B_k = B/trace(B); % Normalize it based on its trace. 
    B_k = B;

    lambda = diag(Lambda);

    [Vx,lambda_x] = eig(L*Gamma*L');
    [Vt,lambda_t] = eig(B_k*eye(T)); 
    lambda_x = spdiags(diag(lambda_x),0,length(lambda_x),length(lambda_x));
    lambda_t = spdiags(diag(lambda_t),0,length(lambda_t),length(lambda_t));

    for  Eff_index_i = 1:M
        for Eff_index_j =1:T
            temp = ( (lambda_x(Eff_index_i,Eff_index_i) * lambda_t(Eff_index_j,Eff_index_j)) + lambda(Eff_index_i)); 
            P(Eff_index_j,Eff_index_i) = 1/ temp;              
        end
    end
    for g = 1:G
        mu_eff(((g-1)*block_length)+1:g*block_length,:)= Vt * lambda_t* (P .* (Vt'* Y_block(:,:,g) * Vx)) * Vx' * L * Gamma';
        X_trial(:,:,g) = Vt * lambda_t* (P .* (Vt'* Y_block(:,:,g) * Vx)) * Vx' * L * Gamma';
    end
    X_total =  mu_eff';             

    % ================ Calculate the Log-Likelihood ==================	
	logdet_SigmaY = logdet(real(SigmaY_estimated));
	assert(isreal(logdet_SigmaY),'cov(Y) is not positive !');
    energy = (logdet_SigmaY) + abs(trace(C_y*SigmaY_inv)); % log likelihood of the Gaussian density
    Energy(itr_count) = energy;  
    
	if print_figures 
        figure(2);
        plot((1:itr_count),Energy(1:itr_count));
        title(['Neg-LogLikelihood: ' int2str(itr_count) ' / ' int2str(Max_num_iterations)]);
        xlabel('iteration');
        set(gca(),'XLim',[0 itr_count]);
        drawnow
    end
 
 	itr_count = itr_count + 1; 
    
    if (size(X_total) == size(X_old))
        Converge_X = norm(X_old-X_total,'fro') < EPSILON;
        dX_new = max(max(abs(X_old - X_total)));   
        error_X = norm(X_total-X_old,'fro')^2/norm(X_old,'fro')^2;

        if (PRINT==1) 
            fprintf('dX: %d \n \n',norm(X_old-X_total,'fro'));
            fprintf('dX_new: %d \n',dX_new);
            disp(['X change: ',num2str(max(max(abs(error_X))))]);
        end
        
        if dX_new < EPSILON || (itr_count > Max_num_iterations)
            Converge_X = true; 
            X_itr_est_error(itr_count) = error_X; 
        end
    end
 
	X_old = X_total; 

	eps1=1e-8;
    plot_on = 1; 
    iem = itr_count; 
    nem = Max_num_iterations; 
    like(iem)=-.5*(sum(log(max(d,eps1)))+M*log(2*pi))-.5*sum(sum(C_y.*SigmaY_inv));
%     if(plot_on)
%         figure(3)
%         subplot(2,2,1);plot((1:iem),like(1:iem));
%         title(['Likelihood: ' int2str(iem) ' / ' int2str(nem)]);
%         xlabel('iteration');
%         set(gca(),'XLim',[0 iem]);
%         subplot(2,2,2);
%         gamma_est = zeros(ND,1);
%         gamma_est(keep_list,1) = p_source; 
%         amp_est = sum(reshape(gamma_est,D,N),1);   %voxel power
% 		plot((1:N),amp_est,'r');
%         title(['Voxel power: ' num2str(N) ' / ' num2str(N)]);
%         xlabel('voxel index');
%         set(gca(),'XLim',[1 N]);
%         drawnow
%     end
    
	if print_figures 
	    figure(1)
        gamma_est = zeros(ND,1);
        gamma_est(keep_list,1) = p_source; 
        amp_est = sum(reshape(gamma_est,D,N),1);   %voxel power

%       amp_est = sum(reshape(p_source,D,N),1);   %voxel power

% % 		amp_est = sqrt(sum(real(X_total).^2, 2));
% 		amp_est = amp_est ./ sum(amp_est);
% 		plot((1:N),amp_est,'r');
% 		xlabel('voxel index');
% 		set(gca(),'XLim',[1 N]);
% 		title('\fontsize{16}Estimation') 
% 		drawnow

		
% 		figure(2)
% 		subplot(1,2,2);
% 		imagesc(real(Lambda))
% 		title('\fontsize{16}Estimation') 
% 		axis('equal')
% 		axis('tight')
% 		colorbar
% 		drawnow
%         
%         figure(3)
% 		subplot(1,2,2);
% 		imagesc(real(Gamma))
% 		title('\fontsize{16}Estimation') 
% 		axis('equal')
% 		axis('tight')
% 		colorbar
% 		drawnow
        
	end
end

% Expand hyperparameters 
% gamma_est = real(diag(Gamma));  
% Lambda =  real(Lambda);

gamma_ind = sort(keep_list);
gamma_est = zeros(ND,1);
gamma_est(keep_list,1) = p_source;  

% gamma_est = diag(Gamma);  

SigmaY_estimated = (L * Gamma * L') + (Lambda);
W = Gamma * L' * inv(SigmaY_estimated);
X = X_total;


if (PRINT) 
    fprintf('\nFinish running ...\n'); 
end



return;

end

