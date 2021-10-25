clc
clear
close all

% Setting directoris. 

current_dir = pwd;
cd(current_dir);
addpath('..\..\utils\export_fig\')
addpath('..\..\utils\')

load('..\..\data\data1D.mat', 'L', 'D', 'loc')	

% stdnoise = 0.01;
stdnoise = 0.001;

plot_convergence = true;
% Different number of Blocks or number of trials. 
% Ngrid = 10:10:50;
Ngrid = floor(linspace(5,200,5));
% Ngrid = 100;

repts = 2;
itrN = 0;

for G = Ngrid
    fprintf('Number of temporal blocks or trials:%d\n',G);
    itrN = itrN + 1;
    for i = 1:repts

        % Setting the ground truth (GT) parameters for loading the lead
        % fieled matrix
        N = 2004;   % Number of voxels in source space 
        T = 10;     % Number of time samples in each block
        M = 58;     % Number of sensors

        load('..\..\data\data1D.mat', 'L', 'D', 'loc')	
%         load('data1D.mat', 'L', 'D', 'loc')	

%         % Random Gaussian matrix
%         M = 20;
%         N = 200;
%         T = 10;
%         L = randn(M,N); 
    
        % Generate the temporal correlation matrix with Toeplitz structure
        % B_0 = cov_mat_gen(T,'toeplitz'); 

        % Full-structural 
        % [B_0_root, B_0]  = randpsd(T); 
        [B_0_root, B_0]  = randpsd(T,T,0.5); 


        % Four different modes for generating the spatial correlation matrix
        % gamma_mode = 'identity'; 
        % gamma_mode = 'sparse'; 
        % gamma_mode = 'random';
        % gamma_mode = 'diagonal';

        [Gamma_0,indice] = cov_mat_gen(N,'sparse', 1); 
        % [Gamma_0] = cov_mat_gen(N,'diagonal', 1); 
%         [Gamma_0_root,Gamma_0]  = randpsd(N,N,0.5);

        % Diagonal
        Lambda_0 = cov_mat_gen(M,'diagonal', stdnoise^2);
        
        % Geodesic 
%       [Lambda_0_root, Lamda_0]  = randpsd(T,T,0.5);

        % temporal correlation matrix for noise
        [Epsilon_0_root, Epsilon_0]  = randpsd(T,T,0.5); 
        Epsilon_0 = B_0;
        
        % Covarinace matrix in source space: Kronocker product of spatial and
        % temporal correlation matrices. 
        Sigma_0 = kron(Gamma_0,B_0);
        
%         Sigma_e = kron(Lambda_0,B_0);
        Sigma_e = kron(Lambda_0,Epsilon_0);

        % Covarinace matrix in source space: Kronocker product of spatial and
        % temporal correlation matrices. 
        Ds = diag(eig(Gamma_0));
        Vs = eye(N,N);    
        [Vt Dt] = eig(B_0);
        Mtime_source = Vt*sqrt(Dt);
        Mspace_source = sqrt(Ds)'*Vs';

        % [Vs Ds] = eig((stdnoise^2 * eye(M)) + (L * Gamma_0 * L')); 
        [Vs Ds] = eig(Lambda_0 + (L * Gamma_0 * L')); 
        [Vt_noise Dt_noise] = eig(Epsilon_0);
        
        Mtime_noise = Vt_noise*sqrt(Dt_noise);
        Mspace_sensor = sqrt(Ds)'*Vs';		

        % Covarinace matrix in sensor space: Sigma_Y_temp = Kronocker product of Sigma_y and
        % temporal correlation matrices. Please check Eq. (14) of the T-BSI paper. 

        rng(1);
        rng('default');

        % % == using approximation == %
        SigmaY = kron(Lambda_0 + (L * Gamma_0 * L'), B_0);
        % SigmaY = kron((stdnoise^2 * eye(M)) + (L * Gamma_0 * L'), B_0);
        % SigmaY_test = kron((stdnoise^2 * eye(M)) + (L * Gamma_0 * L'), B_0);

        % % == without approximation using X == % 
        % I_T= spdiags(ones(T),0,T,T);
        % I_MT= spdiags(ones(M*T),0,M*T,M*T);
        % 
        % L_prime = kron(L,I_T);
        % LSigma0L = zeros(M*T);
        % 
        % for i = 1 : N
        %     Sigma0 =  Sigma_0((i-1)*T+1: i*T,(i-1)*T+1: i*T); 
        %     LSigma0L = LSigma0L + L_prime(:, (i-1)*T+1: i*T ) * Sigma0 * L_prime(:, (i-1)*T+1: i*T )';
        % end
        % SigmaY = LSigma0L + (stdnoise^2 * I_MT);


        % % Check the error of the approximation
        % error_SigmaY = norm(SigmaY_test/trace(SigmaY_test)-SigmaY/trace(SigmaY),'fro')^2/norm(SigmaY /trace(SigmaY),'fro')^2;
        % fprintf('Approximation error in calcuating covariance matrix in the sensor space covariance:%1.5f', error_SigmaY) 

%         rng(1);
%         % generate samples from a normal distribution 
%         X_total = mvnrnd(zeros(N*T,1),Sigma_0,G)'; 
%         % generate samples Xg
% 
%         % Source Space 
%         for g = 1:G
%             X(:,:,g) = reshape(X_total(:,g),T,N);
%             Y(:,:,g) = L * X(:,:,g)'; 
%         end
%         Y_total_gen = reshape(Y,M*T,G);
%         Y_total_gen_noisy =  Y_total_gen + ((stdnoise^2) * randn(size(Y_total_gen))); 
%         Y_total_sensor = Y_total_gen_noisy; 

        % Sensor Space 
        rng(1);
		rng('default');
%         X_0 = mvnrnd(zeros(N*T,1),Sigma_0,G)'; 
%         Y_total = mvnrnd(zeros(M*T,1),SigmaY,G)'; 
        
        E_total = mvnrnd(zeros(M*T,1),Sigma_e,G)';
        clear z_emp1 z_emp3
        for isamp = 1:G
          z_emp1(:, isamp) = vec(Mtime_source*randn(T, N)*Mspace_source);
          z_emp2(:, isamp) = vec(Mtime_noise*randn(T, M)*Mspace_sensor);
        end

        X_0 = z_emp1; 
         
        for g = 1:G
            X(:,:,g) = reshape(X_0(:,g),T,N);

%             %% Adding the noise to each trial 
%             Y_trial_vec = L * X(:,:,g)'; 
%             norm_signal = norm(Y_trial_vec, 'fro');
%             E_trial = mvnrnd(zeros(M*T,1),Sigma_e,1)'; 
%             noise = reshape(E_trial,M,T);
%             norm_noise = norm(noise, 'fro'); 
%             noise = noise ./ norm_noise; 
% 
%             % SNR (dB) based on the energy ration between signal and noise
%             alpha = 0.9; 
%             SNR_value = 20*log10(alpha/(1-alpha));
%             Y_trial(:,:,g) = Y_trial_vec + (1-alpha)*noise*norm_signal/alpha;
%             EEG_baseline(:,:,g) = (1-alpha)*noise*norm_signal/alpha;
% 
%             scale = (1-alpha)*norm_signal/(alpha*norm_noise); 
%             Lambda_0 = (scale^2) * Lambda_0; 
%             Sigma_e = (scale^2) * Sigma_e;
%             stdnoise = scale * stdnoise; 
        end
%         Y_total = reshape(Y_trial,[],G);
        
%         X_avr = mean(X,3)';
%         X_0 = reshape(X_0,N,[]);
        
          X_0 = reshape(permute(X,[2,1,3]),N,[]);
%         Y_total = reshape(permute(Y_trial,[3,1,2]),G,[])';
        
        
%% Adding the noise with respective SNR 
%         Y_nonoise = L * X_0; 
%         norm_signal = norm(Y_nonoise, 'fro');
%         noise = E_total; 
%         norm_noise = norm(noise, 'fro'); 
%         noise = noise ./ norm_noise; 
% 
%         % SNR (dB) based on the energy ration between signal and noise
%         alpha = 0.9; 
%         SNR_value = 20*log10(alpha/(1-alpha));
%         Y_signal = noise + (1-alpha)*noise*norm_signal/alpha;
%         EEG_baseline = (1-alpha)*noise*norm_signal/alpha;
%         EEG_baseline_avr = mean(permute(reshape(permute(EEG_baseline,[2 1]),[G,M,T]),[2,3,1]),3);
%         scale = (1-alpha)*norm_signal/(alpha*norm_noise); 
%         Lambda_0 = (scale^2) * Lambda_0; 
%         Sigma_e = (scale^2) * Sigma_e;
%         stdnoise = scale * stdnoise; 
%         
%         Y_total = Y_signal;   
%% 

        
        Y_total = z_emp2;  
        noise_cov = Lambda_0; 

        %% A simple check        
        %       % Check the difference between source and sensor space       
        %       error_source_sensor = norm(Y_total-Y_total_sensor,'fro')^2/norm(Y_total_sensor,'fro')^2;
        %       fprintf('Error in generatng sensor data using X and Y:%3.2f', error_source_sensor) 
 
        
%         %% estimate noise variance
%         EEG_baseline_avr = mean(EEG_baseline,3);
%         plotflag=1;
%         nl= 2;
%         nem = 100;       
%         [b,lam,sig,yc,cy,bet,weight,mlike,xubar]=nut_reg_vbfa(EEG_baseline_avr,nl,nem,plotflag);
%         noise_cov = double(sig);
		%% ---------------------------------------------
		% Normalizing the leadfiled
		DW = sqrt(sum(L.^2));
		L = L*diag(1./DW);
		
		%% 
        tic
        [X_total_toeplitz,X_trials_toeplitz,Sigma_toeplitz, Gamma_toeplitz, B_toeplitz,Gamma_error_toeplitz, B_error_toeplitz,f_loglikelihood_toeplitz]= ...
            CovarianceEstimator_sensor_SigmaY_X_stable(Y_total,L,SigmaY,noise_cov,'sum_of_rank_one', 'Toeplitz');
        time_toeplitz = toc; 
        
%         tic
%         [X_total_toeplitz,X_trials_toeplitz,Sigma_toeplitz, Gamma_toeplitz, B_toeplitz,Gamma_error_toeplitz, B_error_toeplitz,f_loglikelihood_toeplitz]= ...
%             CovarianceEstimator_sensor_SigmaY_X_stable(Y_total,L,SigmaY,Lambda_0,'Geodesic', 'Toeplitz');
%         time_toeplitz = toc; 
        
        
        X_total_toeplitz = diag(1./DW) * real(X_total_toeplitz);  

        corr_B_toeplitz = corr(real(B_toeplitz(:)),B_0(:));
        NMSE_B_toeplitz = norm(B_0-real(B_toeplitz),'fro')^2/norm(B_0,'fro')^2;
        
        corr_Gamma_toeplitz = corr(real(Gamma_toeplitz(:)),Gamma_0(:));
        NMSE_Gamma_toeplitz = norm(Gamma_0-real(Gamma_toeplitz),'fro')^2/norm(Gamma_0,'fro')^2;

		%%         
        tic
        [X_total_geodesic,X_trials_geodesic,Sigma_geodesic, Gamma_geodesic, B_geodesic,Gamma_error_geodesic, B_error_geodesic,f_loglikelihood_geodesic]= ... 
            CovarianceEstimator_sensor_SigmaY_X_stable(Y_total,L,SigmaY,noise_cov,'sum_of_rank_one', 'Geodesic');
        time_geodesic = toc;
        
%         tic
%         [X_total_geodesic,X_trials_geodesic,Sigma_geodesic, Gamma_geodesic, B_geodesic,Gamma_error_geodesic, B_error_geodesic,f_loglikelihood_geodesic]= ... 
%             CovarianceEstimator_sensor_SigmaY_X_stable(Y_total,L,SigmaY,Lambda_0,'Geodesic', 'Geodesic');
%         time_geodesic = toc; 

        X_total_geodesic = diag(1./DW) * real(X_total_geodesic);  

        corr_B_geodesic = corr(real(B_geodesic(:)),B_0(:));
        NMSE_B_geodesic = norm(B_0-real(B_geodesic),'fro')^2/norm(B_0,'fro')^2;
        
        corr_Gamma_geodesic = corr(real(Gamma_geodesic(:)),Gamma_0(:));
        NMSE_Gamma_geodesic = norm(Gamma_0-real(Gamma_geodesic),'fro')^2/norm(Gamma_0,'fro')^2;


%%      Running Champagne 
        for g = 1:G
            Y(:,:,g)=reshape(Y_total(:,g),T,M);
        end
%         Y_ave = mean(Y,3)';
%         Y_champ = Y_ave; 
        Y_champ = reshape(Y,M,[]); 
        
        tic 
		[gamma,x,w,c,Gamma_error_champ,f_loglikelihood_champ]= awsm_champ(Y_champ,L,noise_cov,150,1,1,1); 
        time_champ = toc;
        X_champ = diag(1./DW) * real(x);  
		Gamma_champ = diag(squeeze(gamma(1,1,:))); 
        Sigma_champ = c; 
        
        
        %% Calculating EMD and time course Error
        [EMD_toeplitz,Corr_toeplitz,EUCL_toeplitz] = perf_emd(X_total_toeplitz, X_0, D, indice);
        [EMD_geodesic,Corr_geodesic,EUCL_geodesic] = perf_emd(X_total_geodesic, X_0, D, indice);
        [EMD_champ,Corr_champ,EUCL_champ] = perf_emd(X_champ, X_0, D, indice);
        

        %% Calculating F1-Score
        F1measure_toeplitz = calc_F1measure(X_total_toeplitz,X_0);
        F1measure_geodesic = calc_F1measure(X_total_geodesic,X_0);
        F1measure_champ = calc_F1measure(X_champ,X_0);
        
        %% Calculating MSE
        MSE_toeplitz = (norm(X_total_toeplitz-X_0,'fro')/norm(X_0,'fro'))^2;
        MSE_geodesic = (norm(X_total_geodesic-X_0,'fro')/norm(X_0,'fro'))^2;
        MSE_champ = (norm(X_champ-X_0,'fro')/norm(X_0,'fro'))^2;

        
        fprintf('Toeplitz: MSE = %3.2f; B-Corr = %4.3f; F1 = %4.3f; EMD = %4.3f; CORR = %4.3f; EUCL = %4.3f; Time = %4.3f \n',...
            MSE_toeplitz, corr_B_toeplitz, F1measure_toeplitz,  EMD_toeplitz,  Corr_toeplitz, EUCL_toeplitz, time_toeplitz);
    
        fprintf('Geodesic: MSE = %3.2f; B-Corr = %4.3f; F1 = %4.3f; EMD = %4.3f; CORR = %4.3f; EUCL = %4.3f; Time = %4.3f \n',...
            MSE_geodesic, corr_B_geodesic, F1measure_geodesic,  EMD_geodesic,  Corr_geodesic, EUCL_geodesic, time_geodesic);
        
        fprintf('Champagne: MSE = %3.2f; F1 = %4.3f; EMD = %4.3f; CORR = %4.3f; EUCL = %4.3f; Time = %4.3f \n',...
            MSE_champ, F1measure_champ,  EMD_champ,  Corr_champ, EUCL_champ, time_champ);

     %%   
%         if plot_convergence
% %%        
%             % figure
%             figure('units','normalized','outerposition',[0 0 0.65 0.85])
%             % plot(Ngrid*T,err_geodesic_structure,'b+-',Ngrid*T,err_toeplitz_structure,'r.-.');
% 
%             h_geodesic = plot(squeeze(B_error_geodesic));
%             LineStyle = '--';
%             LineWidth = 4;    
%             Color = [1 0 1];      
%             Marker =  'd';   
%             MarkerSize = 9;      
%             MarkerEdgeColor = [1 0 1]; 
%             MarkerFaceColor = [1 0 1];
%             set(h_geodesic                                  , ...
%             'LineStyle'       , LineStyle         , ... 
%             'LineWidth'       , LineWidth         , ...
%             'Color'           , Color             , ...
%             'Marker'          , Marker            , ...
%             'MarkerSize'      , MarkerSize    	, ...
%             'MarkerEdgeColor' , MarkerEdgeColor   , ...
%             'MarkerFaceColor' , MarkerFaceColor );
% 
%             hold on
%             h_toeplitz = plot(squeeze(B_error_toeplitz));
%             LineStyle = '--';
%             LineWidth = 4;         
%             Color = [0 0.5 .5];      
%             Marker =  'o';   
%             MarkerSize = 9;      
%             MarkerEdgeColor = [0 0.5 .5]; 
%             MarkerFaceColor = [0 0.5 .5];
%             set(h_toeplitz                                  , ...
%             'LineStyle'       , LineStyle         , ... 
%             'LineWidth'       , LineWidth         , ...
%             'Color'           , Color             , ...
%             'Marker'          , Marker            , ...
%             'MarkerSize'      , MarkerSize    	, ...
%             'MarkerEdgeColor' , MarkerEdgeColor   , ...
%             'MarkerFaceColor' , MarkerFaceColor );
% 
%             hLegend = legend('Geodesic','Toeplitz');  
%             set(hLegend,'interpreter','tex');
%             set(hLegend,'Location','Northeast');
%             set([hLegend, gca]             , ...
%                 'FontSize'   , 30           );
% 
%             set(gca, ...
%               'Box'         , 'off'     , ...
%               'TickDir'     , 'in'     , ...
%               'YMinorTick'  , 'on'      , ...
%               'TickLength'  , [.01 .05] , ...
%               'YGrid'       , 'on'      , ...
%               'XMinorTick'  , 'on'      , ...
%               'YMinorGrid'  , 'on'      , ...
%               'XGrid'       , 'on'      , ...
%               'XMinorTick'  , 'on'      , ...
%               'TickLength'  , [.01 .05] , ...
%               'FontSize'    , 30        , ... 
%               'LineWidth'   , 1         ); 
%             set(gca, 'YScale', 'log')
%             set(gca, 'XScale', 'log')		
%             set(gcf,'color','w');
% 
%             hXLabel = xlabel('Number of iterations');
%             hYLabel = ylabel('Estimation Error (NMSE)');
%             title('Estimation error for temporal covariance') 
%             set(hXLabel,'Interpreter','tex');
%             set(hYLabel,'Interpreter','tex');
%             set(hXLabel,'FontSize',40);
%             set(hYLabel,'FontSize',40);
% 
%             res = 100; 
%             file_name = ['temporal_NMSE_Gamma']; 
%             savefig(file_name); 
%             export_fig(file_name, ['-r' num2str(res)], '-a2', '-transparent', '-eps', '-pdf', '-png'); 
%     %% 
%             % figure
%             figure('units','normalized','outerposition',[0 0 0.65 0.85])
%             % plot(Ngrid*T,err_geodesic_structure,'b+-',Ngrid*T,err_toeplitz_structure,'r.-.');
% 
%             h_toeplitz_SumofRankone = plot(squeeze(Gamma_error_toeplitz));
%             LineStyle = '--';
%             LineWidth = 4;    
%             Color = [1 0 1];      
%             Marker =  'd';   
%             MarkerSize = 9;      
%             MarkerEdgeColor = [1 0 1]; 
%             MarkerFaceColor = [1 0 1];
%             set(h_toeplitz_SumofRankone                                  , ...
%             'LineStyle'       , LineStyle         , ... 
%             'LineWidth'       , LineWidth         , ...
%             'Color'           , Color             , ...
%             'Marker'          , Marker            , ...
%             'MarkerSize'      , MarkerSize    	, ...
%             'MarkerEdgeColor' , MarkerEdgeColor   , ...
%             'MarkerFaceColor' , MarkerFaceColor );
% 
%             hold on
%             h_geodesic_SumofRankone = plot(squeeze(Gamma_error_geodesic));
%             LineStyle = '--';
%             LineWidth = 4;         
%             Color = [0 0.5 .5];      
%             Marker =  'o';   
%             MarkerSize = 9;      
%             MarkerEdgeColor = [0 0.5 .5]; 
%             MarkerFaceColor = [0 0.5 .5];
%             set(h_geodesic_SumofRankone                                  , ...
%             'LineStyle'       , LineStyle         , ... 
%             'LineWidth'       , LineWidth         , ...
%             'Color'           , Color             , ...
%             'Marker'          , Marker            , ...
%             'MarkerSize'      , MarkerSize    	, ...
%             'MarkerEdgeColor' , MarkerEdgeColor   , ...
%             'MarkerFaceColor' , MarkerFaceColor );
% 
%             hold on
%             h_champ = plot(squeeze(Gamma_error_champ));
%             LineStyle = '--';
%             LineWidth = 4;         
%             Color = 'r';      
%             Marker =  'p';   
%             MarkerSize = 5;      
%             MarkerEdgeColor = 'r'; 
%             MarkerFaceColor = 'r';
%             set(h_champ                               , ...
%             'LineStyle'       , LineStyle         , ... 
%             'LineWidth'       , LineWidth         , ...
%             'Color'           , Color             , ...
%             'Marker'          , Marker            , ...
%             'MarkerSize'      , MarkerSize    	, ...
%             'MarkerEdgeColor' , MarkerEdgeColor   , ...
%             'MarkerFaceColor' , MarkerFaceColor );
% 
% 
%             hLegend = legend('Toeplitz','Geodesic','Champ');  
%             set(hLegend,'interpreter','tex');
%             set(hLegend,'Location','Northeast');
%             set([hLegend, gca]             , ...
%                 'FontSize'   , 30           );
% 
%             set(gca, ...
%               'Box'         , 'off'     , ...
%               'TickDir'     , 'in'     , ...
%               'YMinorTick'  , 'on'      , ...
%               'TickLength'  , [.01 .05] , ...
%               'YGrid'       , 'on'      , ...
%               'XMinorTick'  , 'on'      , ...
%               'YMinorGrid'  , 'on'      , ...
%               'XGrid'       , 'on'      , ...
%               'XMinorTick'  , 'on'      , ...
%               'TickLength'  , [.01 .05] , ...
%               'FontSize'    , 30        , ... 
%               'LineWidth'   , 1         ); 
%             set(gca, 'YScale', 'log')
%             set(gca, 'XScale', 'log')		
%             set(gcf,'color','w');
% 
%             hXLabel = xlabel('Number of iterations');
%             hYLabel = ylabel('Estimation Error (NMSE)');
%             title('Estimation error for spatial covariance') 
%             set(hXLabel,'Interpreter','tex');
%             set(hYLabel,'Interpreter','tex');
%             set(hXLabel,'FontSize',40);
%             set(hYLabel,'FontSize',40);
% 
%             res = 100; 
%             file_name = ['spatial_NMSE_Gamma']; 
%             savefig(file_name); 
%             export_fig(file_name, ['-r' num2str(res)], '-a2', '-transparent', '-eps', '-pdf', '-png'); 

%     %% 
%             % figure
%             figure('units','normalized','outerposition',[0 0 0.65 0.85])
%             % plot(Ngrid*T,err_geodesic_structure,'b+-',Ngrid*T,err_toeplitz_structure,'r.-.');
% 
%             h_toeplitz_SumofRankone = plot(squeeze(abs(f_loglikelihood_toeplitz)));
%             LineStyle = '--';
%             LineWidth = 4;    
%             Color = [1 0 1];      
%             Marker =  'd';   
%             MarkerSize = 9;      
%             MarkerEdgeColor = [1 0 1]; 
%             MarkerFaceColor = [1 0 1];
%             set(h_toeplitz_SumofRankone                                  , ...
%             'LineStyle'       , LineStyle         , ... 
%             'LineWidth'       , LineWidth         , ...
%             'Color'           , Color             , ...
%             'Marker'          , Marker            , ...
%             'MarkerSize'      , MarkerSize    	, ...
%             'MarkerEdgeColor' , MarkerEdgeColor   , ...
%             'MarkerFaceColor' , MarkerFaceColor );
% 
%             hold on
%             h_geodesic = plot(squeeze(abs(f_loglikelihood_geodesic)));
%             LineStyle = '--';
%             LineWidth = 4;         
%             Color = [0 0.5 .5];      
%             Marker =  'o';   
%             MarkerSize = 9;      
%             MarkerEdgeColor = [0 0.5 .5]; 
%             MarkerFaceColor = [0 0.5 .5];
%             set(h_geodesic                                  , ...
%             'LineStyle'       , LineStyle         , ... 
%             'LineWidth'       , LineWidth         , ...
%             'Color'           , Color             , ...
%             'Marker'          , Marker            , ...
%             'MarkerSize'      , MarkerSize    	, ...
%             'MarkerEdgeColor' , MarkerEdgeColor   , ...
%             'MarkerFaceColor' , MarkerFaceColor );
% 
%             hold on
%             h_champ = plot(squeeze(f_loglikelihood_champ));
%             LineStyle = '--';
%             LineWidth = 4;         
%             Color = 'r';      
%             Marker =  'p';   
%             MarkerSize = 5;      
%             MarkerEdgeColor = 'r'; 
%             MarkerFaceColor = 'r';
%             set(h_champ                               , ...
%             'LineStyle'       , LineStyle         , ... 
%             'LineWidth'       , LineWidth         , ...
%             'Color'           , Color             , ...
%             'Marker'          , Marker            , ...
%             'MarkerSize'      , MarkerSize    	, ...
%             'MarkerEdgeColor' , MarkerEdgeColor   , ...
%             'MarkerFaceColor' , MarkerFaceColor );
% 
%             hLegend = legend('Toeplitz','Geodesic','Champ');  
%             set(hLegend,'interpreter','tex');
%             set(hLegend,'Location','Northeast');
%             set([hLegend, gca]             , ...
%                 'FontSize'   , 30           );
% 
%             set(gca, ...
%               'Box'         , 'off'     , ...
%               'TickDir'     , 'in'     , ...
%               'YMinorTick'  , 'on'      , ...
%               'TickLength'  , [.01 .05] , ...
%               'YGrid'       , 'on'      , ...
%               'XMinorTick'  , 'on'      , ...
%               'YMinorGrid'  , 'on'      , ...
%               'XGrid'       , 'on'      , ...
%               'XMinorTick'  , 'on'      , ...
%               'TickLength'  , [.01 .05] , ...
%               'FontSize'    , 30        , ... 
%               'LineWidth'   , 1         ); 
%             set(gca, 'YScale', 'log')
%             set(gca, 'XScale', 'log')		
%             set(gcf,'color','w');
% 
%             hXLabel = xlabel('Number of iterations');
%             hYLabel = ylabel('Loss');
%             % title('Estimation error for spatial covariance') 
%             set(hXLabel,'Interpreter','tex');
%             set(hYLabel,'Interpreter','tex');
%             set(hXLabel,'FontSize',40);
%             set(hYLabel,'FontSize',40);
% 
%             res = 100; 
%             file_name = ['NegLogLoss']; 
%             savefig(file_name); 
%             export_fig(file_name, ['-r' num2str(res)], '-a2', '-transparent', '-eps', '-pdf', '-png'); 
%         end
        
        corr_B_structure_toeplitz(itrN,i) = corr(real(B_toeplitz(:)),B_0(:));
        NMSE_B_structure_toeplitz(itrN,i) = norm(B_0-real(B_toeplitz),'fro')^2/norm(B_0,'fro')^2;            
        corr_Gamma_structure_toeplitz(itrN,i) = corr(real(Gamma_toeplitz(:)),Gamma_0(:));
        NMSE_Gamma_structure_toeplitz(itrN,i) = norm(Gamma_0-real(Gamma_toeplitz),'fro')^2/norm(Gamma_0,'fro')^2;
        
        
        MSE_structure_toeplitz(itrN,i) = MSE_toeplitz;
        EMD_structure_toeplitz(itrN,i) = EMD_toeplitz; 
        Corr_structure_toeplitz(itrN,i) = Corr_toeplitz;
        EUCL_structure_toeplitz(itrN,i) = EUCL_toeplitz;
        F1_structure_toeplitz(itrN,i) =  F1measure_toeplitz;
        
        corr_B_structure_geodesic(itrN,i) = corr(real(B_geodesic(:)),B_0(:));
        NMSE_B_structure_geodesic(itrN,i) = norm(B_0-real(B_geodesic),'fro')^2/norm(B_0,'fro')^2;
        corr_Gamma_structure_geodesic(itrN,i) = corr(real(Gamma_geodesic(:)),Gamma_0(:));
        NMSE_Gamma_structure_geodesic(itrN,i) = norm(Gamma_0-real(Gamma_geodesic),'fro')^2/norm(Gamma_0,'fro')^2;

        MSE_structure_geodesic(itrN,i) = MSE_geodesic;
        EMD_structure_geodesic(itrN,i) = EMD_geodesic; 
        Corr_structure_geodesic(itrN,i) = Corr_geodesic;
        EUCL_structure_geodesic(itrN,i) = EUCL_geodesic;
        F1_structure_geodesic(itrN,i) =  F1measure_geodesic;
        
        
        MSE_structure_champ(itrN,i) = MSE_champ;
        EMD_structure_champ(itrN,i) = EMD_champ; 
        Corr_structure_champ(itrN,i) = Corr_champ;
        EUCL_structure_champ(itrN,i) = EUCL_champ;
        F1_structure_champ(itrN,i) =  F1measure_champ;

% % 		err_geodesic_structure(itrN,i) = norm(SigmaY/trace(SigmaY)-Sigma_geodesic/trace(Sigma_geodesic),'fro')^2/norm(SigmaY/trace(SigmaY),'fro')^2;
% 		err_geodesic_structure_Gamma(itrN,i) = norm(Gamma_0/trace(Gamma_0)-Gamma_geodesic/trace(Gamma_geodesic),'fro')^2/norm(Gamma_0/trace(Gamma_0),'fro')^2;
% 		err_geodesic_structure_B(itrN,i) = norm(B_0/trace(B_0)-B_geodesic/trace(B_geodesic),'fro')^2/norm(B_0/trace(B_0),'fro')^2;
% 		
% % 		err_toeplitz_structure(itrN,i) = norm(SigmaY/trace(SigmaY)-Sigma_toeplitz/trace(Sigma_toeplitz),'fro')^2/norm(SigmaY/trace(SigmaY),'fro')^2;
% 		err_toeplitz_structure_Gamma(itrN,i) = norm(Gamma_0/trace(Gamma_0)-Gamma_toeplitz/trace(Gamma_toeplitz),'fro')^2/norm(Gamma_0/trace(Gamma_0),'fro')^2;
% 		err_toeplitz_structure_B(itrN,i) = norm(B_0/trace(B_0)-B_toeplitz/trace(B_toeplitz),'fro')^2/norm(B_0/trace(B_0),'fro')^2;
% 
%         err_champ(itrN,i) = norm(SigmaY/trace(SigmaY)-Sigma_champ/trace(Sigma_champ),'fro')^2/norm(SigmaY/trace(SigmaY),'fro')^2;
% 		err_champ_Gamma(itrN,i) = norm(Gamma_0/trace(Gamma_0)-Gamma_champ/trace(Gamma_champ),'fro')^2/norm(Gamma_0/trace(Gamma_0),'fro')^2;

    end
end
  
 
corr_B_structure_toeplitz = sum(corr_B_structure_toeplitz,2)/repts;
NMSE_B_structure_toeplitz = sum(NMSE_B_structure_toeplitz,2)/repts;
corr_Gamma_structure_toeplitz = sum(corr_Gamma_structure_toeplitz,2)/repts;
NMSE_Gamma_structure_toeplitz = sum(NMSE_Gamma_structure_toeplitz,2)/repts;

MSE_structure_toeplitz = sum(MSE_structure_toeplitz,2)/repts;
EMD_structure_toeplitz = sum(EMD_structure_toeplitz,2)/repts;
Corr_structure_toeplitz = sum(Corr_structure_toeplitz,2)/repts;
EUCL_structure_toeplitz = sum(EUCL_structure_toeplitz,2)/repts;
F1_structure_toeplitz = sum(F1_structure_toeplitz,2)/repts;

corr_B_structure_geodesic = sum(corr_B_structure_geodesic,2)/repts;
NMSE_B_structure_geodesic = sum(NMSE_B_structure_geodesic,2)/repts;
corr_Gamma_structure_geodesic = sum(corr_Gamma_structure_geodesic,2)/repts;
NMSE_Gamma_structure_geodesic = sum(NMSE_Gamma_structure_geodesic,2)/repts;

MSE_structure_geodesic = sum(MSE_structure_geodesic,2)/repts;
EMD_structure_geodesic = sum(EMD_structure_geodesic,2)/repts;
Corr_structure_geodesic = sum(Corr_structure_geodesic,2)/repts;
EUCL_structure_geodesic = sum(EUCL_structure_geodesic,2)/repts;
F1_structure_geodesic = sum(F1_structure_geodesic,2)/repts;


MSE_structure_champ = sum(MSE_structure_champ,2)/repts;
EMD_structure_champ = sum(EMD_structure_champ,2)/repts;
Corr_structure_champ = sum(Corr_structure_champ,2)/repts;
EUCL_structure_champ = sum(EUCL_structure_champ,2)/repts;
F1_structure_champ = sum(F1_structure_champ,2)/repts;


% err_geodesic_structure_Gamma = sum(err_geodesic_structure_Gamma,2)/repts;
% err_toeplitz_structure_Gamma = sum(err_toeplitz_structure_Gamma,2)/repts;
% 
% 
% err_geodesic_structure_B = sum(err_geodesic_structure_B,2)/repts;
% err_toeplitz_structure_B = sum(err_toeplitz_structure_B,2)/repts;
 
%% 
% Ploting the MSE error for comparing the perfomance of two methods for 
% estimating the spatial and temporal covariance matrices

% figure
figure('units','normalized','outerposition',[0 0 0.65 0.85])
% plot(Ngrid*T,err_geodesic_structure,'b+-',Ngrid*T,err_toeplitz_structure,'r.-.');

h_toeplitz_SumofRankone = plot(Ngrid,corr_B_structure_toeplitz);
LineStyle = '--';
LineWidth = 4;    
Color = [1 0 1];      
Marker =  'd';   
MarkerSize = 9;      
MarkerEdgeColor = [1 0 1]; 
MarkerFaceColor = [1 0 1];
set(h_toeplitz_SumofRankone                                  , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , MarkerEdgeColor   , ...
'MarkerFaceColor' , MarkerFaceColor );

hold on
h_geodesic = plot(Ngrid,corr_B_structure_geodesic);
LineStyle = '--';
LineWidth = 4;         
Color = [0 0.5 .5];      
Marker =  'o';   
MarkerSize = 9;      
MarkerEdgeColor = [0 0.5 .5]; 
MarkerFaceColor = [0 0.5 .5];
set(h_geodesic                                  , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , MarkerEdgeColor   , ...
'MarkerFaceColor' , MarkerFaceColor );

hLegend = legend('Toeplitz','Geodesic');  
set(hLegend,'interpreter','tex');
set(hLegend,'Location','Northeast');
set([hLegend, gca]             , ...
    'FontSize'   , 30           );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'YMinorTick'  , 'on'      , ...
  'TickLength'  , [.01 .05] , ...
  'YGrid'       , 'on'      , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorGrid'  , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'XMinorTick'  , 'on'      , ...
  'TickLength'  , [.01 .05] , ...
  'FontSize'    , 30        , ... 
  'LineWidth'   , 1         ); 
set(gcf,'color','w');

hXLabel = xlabel('Number of Trials');
hYLabel = ylabel('Corr');
% title('Estimation error for temporal Cov matrix') 
set(hXLabel,'Interpreter','tex');
set(hYLabel,'Interpreter','tex');
set(hXLabel,'FontSize',40);
set(hYLabel,'FontSize',40);

res = 100; 
file_name = ['spatiotemporal_NMSE_senosr_space']; 
savefig(file_name); 
export_fig(file_name, ['-r' num2str(res)], '-a2', '-transparent', '-eps', '-pdf'); 
res = 150; 
export_fig(file_name,['-r' num2str(res)], '-png'); 
%% 
% Ploting the MSE error for comparing the perfomance of two methods for 
% estimating the spatial covariance matrice 

% figure
figure('units','normalized','outerposition',[0 0 0.65 0.85])
% plot(Ngrid*T,err_geodesic_structure,'b+-',Ngrid*T,err_toeplitz_structure,'r.-.');

h_toeplitz_SumofRankone = plot(Ngrid,NMSE_B_structure_toeplitz);
LineStyle = '--';
LineWidth = 4;    
Color = [1 0 1];      
Marker =  'd';   
MarkerSize = 9;      
MarkerEdgeColor = [1 0 1]; 
MarkerFaceColor = [1 0 1];
set(h_toeplitz_SumofRankone           , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , MarkerEdgeColor   , ...
'MarkerFaceColor' , MarkerFaceColor );

hold on
h_geodesic = plot(Ngrid,NMSE_Gamma_structure_toeplitz);
LineStyle = '--';
LineWidth = 4;         
Color = [0 0.5 .5];      
Marker =  'o';   
MarkerSize = 9;      
MarkerEdgeColor = [0 0.5 .5]; 
MarkerFaceColor = [0 0.5 .5];
set(h_geodesic                   	  , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , MarkerEdgeColor   , ...
'MarkerFaceColor' , MarkerFaceColor );

hLegend = legend('Toeplitz','Geodesic');  
set(hLegend,'interpreter','tex');
set(hLegend,'Location','Northeast');
set([hLegend, gca]             , ...
    'FontSize'   , 30           );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'YMinorTick'  , 'on'      , ...
  'TickLength'  , [.01 .05] , ...
  'YGrid'       , 'on'      , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorGrid'  , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'XMinorTick'  , 'on'      , ...
  'TickLength'  , [.01 .05] , ...
  'FontSize'    , 30        , ... 
  'LineWidth'   , 1         ); 
set(gcf,'color','w');

hXLabel = xlabel('Number of Trials');
hYLabel = ylabel('NMSE');
% title('Estimation error for temporal covarinace') 
set(hXLabel,'Interpreter','tex');
set(hYLabel,'Interpreter','tex');
set(hXLabel,'FontSize',40);
set(hYLabel,'FontSize',40);

res = 100; 
file_name = ['spatiotemporal_NMSE_senosr_space_B']; 
savefig(file_name); 
export_fig(file_name, ['-r' num2str(res)], '-a2', '-transparent', '-eps', '-pdf'); 
res = 150; 
export_fig(file_name,['-r' num2str(res)], '-png'); 


%% EMD
% figure
figure('units','normalized','outerposition',[0 0 0.65 0.85])
% plot(Ngrid*T,err_geodesic_structure,'b+-',Ngrid*T,err_toeplitz_structure,'r.-.');

h_toeplitz_SumofRankone = plot(Ngrid,EMD_toeplitz);
LineStyle = '--';
LineWidth = 4;    
Color = [1 0 1];      
Marker =  'd';   
MarkerSize = 9;      
MarkerEdgeColor = [1 0 1]; 
MarkerFaceColor = [1 0 1];
set(h_toeplitz_SumofRankone           , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , MarkerEdgeColor   , ...
'MarkerFaceColor' , MarkerFaceColor );

hold on
h_geodesic = plot(Ngrid,EMD_geodesic);
LineStyle = '--';
LineWidth = 4;         
Color = [0 0.5 .5];      
Marker =  'o';   
MarkerSize = 9;      
MarkerEdgeColor = [0 0.5 .5]; 
MarkerFaceColor = [0 0.5 .5];
set(h_geodesic                   	  , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , MarkerEdgeColor   , ...
'MarkerFaceColor' , MarkerFaceColor );

hold on
h_champ = plot(Ngrid,EMD_champ);
LineStyle = '--';
LineWidth = 4;         
Color = 'b';      
Marker =  'o';   
MarkerSize = 9;      
MarkerEdgeColor = [0 0.5 .5]; 
MarkerFaceColor = [0 0.5 .5];
set(h_champ                  	  , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , Color   , ...
'MarkerFaceColor' , Color);

hLegend = legend('Toeplitz','Geodesic','Champ');  
set(hLegend,'interpreter','tex');
set(hLegend,'Location','Northeast');
set([hLegend, gca]             , ...
    'FontSize'   , 30           );

% hLegend = legend([h_toeplitz_SumofRankone ;h_geodesic; h_champ], {'Geodesic','Toeplitz','Champagne'});
% hLegend = legend(h_toeplitz_SumofRankone,'Toeplitz',h_geodesic,'Geodesic',h_champ,'Champagne');  
% hLegend = legend( ...
% [h_toeplitz_SumofRankone, ...  
% h_geodesic, ... ;
% h_champ], ...
% 'Toeplitz', ...
% 'Geodesic', ...
% 'Champagne'); 

% set(hLegend,'interpreter','tex');
% set(hLegend,'Location','Northeast');
% set([hLegend, gca]             , ...
%     'FontSize'   , 30           );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'YMinorTick'  , 'on'      , ...
  'TickLength'  , [.01 .05] , ...
  'YGrid'       , 'on'      , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorGrid'  , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'XMinorTick'  , 'on'      , ...
  'TickLength'  , [.01 .05] , ...
  'FontSize'    , 30        , ... 
  'LineWidth'   , 1         ); 
set(gcf,'color','w');

hXLabel = xlabel('Number of Trials');
hYLabel = ylabel('EMD');
% title('Estimation error for temporal covarinace') 
set(hXLabel,'Interpreter','tex');
set(hYLabel,'Interpreter','tex');
set(hXLabel,'FontSize',40);
set(hYLabel,'FontSize',40);

res = 100; 
file_name = ['EMD']; 
savefig(file_name); 
export_fig(file_name, ['-r' num2str(res)], '-a2', '-transparent', '-eps', '-pdf'); 
res = 150; 
export_fig(file_name,['-r' num2str(res)], '-png');


%% Corr
% figure
figure('units','normalized','outerposition',[0 0 0.65 0.85])
% plot(Ngrid*T,err_geodesic_structure,'b+-',Ngrid*T,err_toeplitz_structure,'r.-.');

h_toeplitz_SumofRankone = plot(Ngrid,Corr_toeplitz);
LineStyle = '--';
LineWidth = 4;    
Color = [1 0 1];      
Marker =  'd';   
MarkerSize = 9;      
MarkerEdgeColor = [1 0 1]; 
MarkerFaceColor = [1 0 1];
set(h_toeplitz_SumofRankone           , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , MarkerEdgeColor   , ...
'MarkerFaceColor' , MarkerFaceColor );

hold on
h_geodesic = plot(Ngrid,Corr_geodesic);
LineStyle = '--';
LineWidth = 4;         
Color = [0 0.5 .5];      
Marker =  'o';   
MarkerSize = 9;      
MarkerEdgeColor = [0 0.5 .5]; 
MarkerFaceColor = [0 0.5 .5];
set(h_geodesic                   	  , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , MarkerEdgeColor   , ...
'MarkerFaceColor' , MarkerFaceColor );

hold on
h_champ = plot(Ngrid,Corr_champ);
LineStyle = '--';
LineWidth = 4;         
Color = 'b';      
Marker =  'o';   
MarkerSize = 9;      
MarkerEdgeColor = [0 0.5 .5]; 
MarkerFaceColor = [0 0.5 .5];
set(h_champ                  	  , ...
'LineStyle'       , LineStyle         , ... 
'LineWidth'       , LineWidth         , ...
'Color'           , Color             , ...
'Marker'          , Marker            , ...
'MarkerSize'      , MarkerSize    	, ...
'MarkerEdgeColor' , Color   , ...
'MarkerFaceColor' , Color);

hLegend = legend('Toeplitz','Geodesic','Champ');  
set(hLegend,'interpreter','tex');
set(hLegend,'Location','Northeast');
set([hLegend, gca]             , ...
    'FontSize'   , 30           );

% hLegend = legend([h_toeplitz_SumofRankone ;h_geodesic; h_champ], {'Geodesic','Toeplitz','Champagne'});
% hLegend = legend(h_toeplitz_SumofRankone,'Toeplitz',h_geodesic,'Geodesic',h_champ,'Champagne');  
% hLegend = legend( ...
% [h_toeplitz_SumofRankone, ...  
% h_geodesic, ... ;
% h_champ], ...
% 'Toeplitz', ...
% 'Geodesic', ...
% 'Champagne'); 

% set(hLegend,'interpreter','tex');
% set(hLegend,'Location','Northeast');
% set([hLegend, gca]             , ...
%     'FontSize'   , 30           );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'YMinorTick'  , 'on'      , ...
  'TickLength'  , [.01 .05] , ...
  'YGrid'       , 'on'      , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorGrid'  , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'XMinorTick'  , 'on'      , ...
  'TickLength'  , [.01 .05] , ...
  'FontSize'    , 30        , ... 
  'LineWidth'   , 1         ); 
set(gcf,'color','w');

hXLabel = xlabel('Number of Trials');
hYLabel = ylabel('Time Course Error');
% title('Estimation error for temporal covarinace') 
set(hXLabel,'Interpreter','tex');
set(hYLabel,'Interpreter','tex');
set(hXLabel,'FontSize',40);
set(hYLabel,'FontSize',40);

res = 100; 
file_name = ['time']; 
savefig(file_name); 
export_fig(file_name, ['-r' num2str(res)], '-a2', '-transparent', '-eps', '-pdf'); 
res = 150; 
export_fig(file_name,['-r' num2str(res)], '-png');