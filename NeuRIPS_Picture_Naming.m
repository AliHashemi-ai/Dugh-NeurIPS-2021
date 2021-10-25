clc
clear
close all

% Setting directoris.

current_dir = pwd;
cd(current_dir);
addpath('../../utils/export_fig/')
addpath('../../utils/')

%% Add Paths (nutmeg & spm)
% addpath('/netopt/share/bin/local/bil/nutmeg4.1') % Chang's version
addpath('/netopt/share/bin/local/bil/nutmeg4b1') %Kamalini's suggested ver.
addpath('/netopt/share/bin/local/bil/spm8_nutmeg_lhinkley')

nuts = load('R2084B.mat'); % picture-naming 1st subject
%nuts = load('C3491A_AEF_02'); % AEF 1st subject
% nuts = load('C3504A_AEF_02'); % AEF 2nd subject
rt = nuts.meg.srate;
data = nuts.meg.data;

[sp nc nt] = size(data);
for i=1:nc
    for j =1:nt
        data(:,i,j) = detrend(data(:,i,j));
    end
end

data = filter_cc(data,'firls','bp','auto',1,70,rt,1);


%% set time point and preprocess leadfields
nsp = size(nuts.meg.latency,1);
pre_start=1;
pre_end = find(0==nuts.meg.latency);%-0.05*1200;
% for AEF task
% post_start = pre_end+.05*1200;%+0.00*1200; % Default setting from Chang
% post_end = pre_end+.15*1200;%nsp; %post_start+.1*1200;%nsp;
% for picture-naming
post_start = pre_end+.025*1200;%+0.00*1200;
post_end = pre_end+.150*1200;%nsp; %post_start+.1*1200;%nsp;
% post_end = pre_end+.275*1200;%nsp; %post_start+.1*1200;%nsp;
% for picture-naming data, task lies in (25-275ms)


[nc,nd,nv] = size(nuts.Lp);

LF = reshape(nuts.Lp,nc,nd*nv);
for i=1:size(LF,2)
    LF(:,i) = LF(:,i)./sqrt(sum(LF(:,i).^2));
end
% LF_eLORETA = reshape(LF,nc,nv,nd);
LF_eLORETA_tmp = reshape(LF,nc,nd,nv);
LF_eLORETA = permute(LF_eLORETA_tmp,[1 3 2]);
LF2 = LF;

voxels = nuts.voxels;
voxelsize=nuts.voxelsize;
coreg=nuts.coreg;
bands=[1 70];
srate=nuts.meg.srate;
timepts = nuts.meg.latency(post_start:post_end,1);


% %%
% load('data/data1D.mat', 'L', 'D', 'loc')
% 
% % stdnoise = 0.01;
% stdnoise = 0.001;
% 
% plot_convergence = true;
% % Different number of Blocks or number of trials.
% % Ngrid = 10:10:50;
% Ngrid = floor(linspace(5,200,5));
% % Ngrid = 100;
% 
% repts = 2;
% itrN = 0;
% 
% %     for G = Ngrid
% %     fprintf('Number of temporal blocks or trials:%d\n',G);
% %     itrN = itrN + 1;
% %     for i = 1:repts
% 
% % Setting the ground truth (GT) parameters for loading the lead
% % fieled matrix
% N = 2004;   % Number of voxels in source space
% T = 10;     % Number of time samples in each block
% M = 58;     % Number of sensors
% 
% load('data/data1D.mat', 'L', 'D', 'loc')
% %         load('data1D.mat', 'L', 'D', 'loc')
% 
% %         % Random Gaussian matrix
% %         M = 20;
% %         N = 200;
% %         T = 10;
% %         L = randn(M,N);
% 
% % Generate the temporal correlation matrix with Toeplitz structure
% % B_0 = cov_mat_gen(T,'toeplitz');
% 
% % Full-structural
% % [B_0_root, B_0]  = randpsd(T);
% [B_0_root, B_0]  = randpsd(T,T,0.5);
% 
% 
% % Four different modes for generating the spatial correlation matrix
% % gamma_mode = 'identity';
% % gamma_mode = 'sparse';
% % gamma_mode = 'random';
% % gamma_mode = 'diagonal';
% 
% [Gamma_0,indice] = cov_mat_gen(N,'sparse', 1);
% % [Gamma_0] = cov_mat_gen(N,'diagonal', 1);
% %         [Gamma_0_root,Gamma_0]  = randpsd(N,N,0.5);
% 
% % Diagonal
% Lambda_0 = cov_mat_gen(M,'diagonal', stdnoise^2);
% 
% % Geodesic
% %       [Lambda_0_root, Lamda_0]  = randpsd(T,T,0.5);
% 
% % temporal correlation matrix for noise
% [Epsilon_0_root, Epsilon_0]  = randpsd(T,T,0.5);
% Epsilon_0 = B_0;
% 
% % Covarinace matrix in source space: Kronocker product of spatial and
% % temporal correlation matrices.
% Sigma_0 = kron(Gamma_0,B_0);
% 
% %         Sigma_e = kron(Lambda_0,B_0);
% Sigma_e = kron(Lambda_0,Epsilon_0);
% 
% % Covarinace matrix in source space: Kronocker product of spatial and
% % temporal correlation matrices.
% Ds = diag(eig(Gamma_0));
% Vs = eye(N,N);
% [Vt Dt] = eig(B_0);
% Mtime_source = Vt*sqrt(Dt);
% Mspace_source = sqrt(Ds)'*Vs';
% 
% % [Vs Ds] = eig((stdnoise^2 * eye(M)) + (L * Gamma_0 * L'));
% [Vs Ds] = eig(Lambda_0 + (L * Gamma_0 * L'));
% [Vt_noise Dt_noise] = eig(Epsilon_0);
% 
% Mtime_noise = Vt_noise*sqrt(Dt_noise);
% Mspace_sensor = sqrt(Ds)'*Vs';
% 
% % Covarinace matrix in sensor space: Sigma_Y_temp = Kronocker product of Sigma_y and
% % temporal correlation matrices. Please check Eq. (14) of the T-BSI paper.
% 
% rng(1);
% rng('default');
% 
% % % == using approximation == %
% SigmaY = kron(Lambda_0 + (L * Gamma_0 * L'), B_0);
% % SigmaY = kron((stdnoise^2 * eye(M)) + (L * Gamma_0 * L'), B_0);
% % SigmaY_test = kron((stdnoise^2 * eye(M)) + (L * Gamma_0 * L'), B_0);
% 
% % % == without approximation using X == %
% % I_T= spdiags(ones(T),0,T,T);
% % I_MT= spdiags(ones(M*T),0,M*T,M*T);
% %
% % L_prime = kron(L,I_T);
% % LSigma0L = zeros(M*T);
% %
% % for i = 1 : N
% %     Sigma0 =  Sigma_0((i-1)*T+1: i*T,(i-1)*T+1: i*T);
% %     LSigma0L = LSigma0L + L_prime(:, (i-1)*T+1: i*T ) * Sigma0 * L_prime(:, (i-1)*T+1: i*T )';
% % end
% % SigmaY = LSigma0L + (stdnoise^2 * I_MT);
% 
% 
% % % Check the error of the approximation
% % error_SigmaY = norm(SigmaY_test/trace(SigmaY_test)-SigmaY/trace(SigmaY),'fro')^2/norm(SigmaY /trace(SigmaY),'fro')^2;
% % fprintf('Approximation error in calcuating covariance matrix in the sensor space covariance:%1.5f', error_SigmaY)
% 
% %         rng(1);
% %         % generate samples from a normal distribution
% %         X_total = mvnrnd(zeros(N*T,1),Sigma_0,G)';
% %         % generate samples Xg
% %
% %         % Source Space
% %         for g = 1:G
% %             X(:,:,g) = reshape(X_total(:,g),T,N);
% %             Y(:,:,g) = L * X(:,:,g)';
% %         end
% %         Y_total_gen = reshape(Y,M*T,G);
% %         Y_total_gen_noisy =  Y_total_gen + ((stdnoise^2) * randn(size(Y_total_gen)));
% %         Y_total_sensor = Y_total_gen_noisy;
% 
% % Sensor Space
% rng(1);
% rng('default');
% %         X_0 = mvnrnd(zeros(N*T,1),Sigma_0,G)';
% %         Y_total = mvnrnd(zeros(M*T,1),SigmaY,G)';
% 
% E_total = mvnrnd(zeros(M*T,1),Sigma_e,G)';
% clear z_emp1 z_emp3
% for isamp = 1:G
%     %           z_emp1(:, isamp) = vec(Mtime_source*randn(T, N)*Mspace_source);
%     %           z_emp2(:, isamp) = vec(Mtime_noise*randn(T, M)*Mspace_sensor);
%     z_emp1(:, isamp) = reshape((Mtime_source*randn(T, N)*Mspace_source)',[],1);
%     z_emp2(:, isamp) = reshape((Mtime_noise*randn(T, M)*Mspace_sensor)',[],1);
% end
% 
% X_0 = z_emp1;
% 
% for g = 1:G
%     X(:,:,g) = reshape(X_0(:,g),T,N);
%     
%     %             %% Adding the noise to each trial
%     %             Y_trial_vec = L * X(:,:,g)';
%     %             norm_signal = norm(Y_trial_vec, 'fro');
%     %             E_trial = mvnrnd(zeros(M*T,1),Sigma_e,1)';
%     %             noise = reshape(E_trial,M,T);
%     %             norm_noise = norm(noise, 'fro');
%     %             noise = noise ./ norm_noise;
%     %
%     %             % SNR (dB) based on the energy ration between signal and noise
%     %             alpha = 0.9;
%     %             SNR_value = 20*log10(alpha/(1-alpha));
%     %             Y_trial(:,:,g) = Y_trial_vec + (1-alpha)*noise*norm_signal/alpha;
%     %             EEG_baseline(:,:,g) = (1-alpha)*noise*norm_signal/alpha;
%     %
%     %             scale = (1-alpha)*norm_signal/(alpha*norm_noise);
%     %             Lambda_0 = (scale^2) * Lambda_0;
%     %             Sigma_e = (scale^2) * Sigma_e;
%     %             stdnoise = scale * stdnoise;
% end
% %         Y_total = reshape(Y_trial,[],G);
% 
% %         X_avr = mean(X,3)';
% %         X_0 = reshape(X_0,N,[]);
% 
% X_0 = reshape(permute(X,[2,1,3]),N,[]);
% %         Y_total = reshape(permute(Y_trial,[3,1,2]),G,[])';
% 
% 
% %% Adding the noise with respective SNR
% %         Y_nonoise = L * X_0;
% %         norm_signal = norm(Y_nonoise, 'fro');
% %         noise = E_total;
% %         norm_noise = norm(noise, 'fro');
% %         noise = noise ./ norm_noise;
% %
% %         % SNR (dB) based on the energy ration between signal and noise
% %         alpha = 0.9;
% %         SNR_value = 20*log10(alpha/(1-alpha));
% %         Y_signal = noise + (1-alpha)*noise*norm_signal/alpha;
% %         EEG_baseline = (1-alpha)*noise*norm_signal/alpha;
% %         EEG_baseline_avr = mean(permute(reshape(permute(EEG_baseline,[2 1]),[G,M,T]),[2,3,1]),3);
% %         scale = (1-alpha)*norm_signal/(alpha*norm_noise);
% %         Lambda_0 = (scale^2) * Lambda_0;
% %         Sigma_e = (scale^2) * Sigma_e;
% %         stdnoise = scale * stdnoise;
% %
% %         Y_total = Y_signal;
% %%
% 
% 
% Y_total = z_emp2;
% noise_cov = Lambda_0;

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

%% Preprocessing & estimate noise variance
num_trials = [40 10 20 60 100];
i_num = 1;
rand_sort = randperm(size(data,3));
%mkdir(nameresult)
randnum = 0;

indx_trial = rand_sort(1:num_trials(i_num));

data_chose = data(:,:,indx_trial);

ypre = permute(data_chose(pre_start:pre_end,:,:),[2 1 3]);
% ypost = permute(data_chose(post_start:post_end,:,:),[2 1 3]);
ypost = data_chose(post_start:post_end,:,:);
ypost_reshape = reshape(ypost,size(ypost,1)*size(ypost,2),size(ypost,3)); 
ypost_test = zeros(size(ypost,1),size(ypost,2),size(ypost,3));
for g = 1:size(ypost,3)
    ypost_test(:,:,g)=reshape(ypost_reshape(:,g),size(ypost,1),size(ypost,2));
end

% Trial concatenation for picture-naming
%         ypre_trial_1=ypre(:,:,1);
%         ypre_trial_2=ypre(:,:,2);

% ypost = reshape(ypost, size(ypost,1),size(ypost,2)*size(ypost,3));
ypre = reshape(ypre, size(ypre,1),size(ypre,2)*size(ypre,3));

%         ypost_reshape_test = reshape(ypost,size(ypost,1),size(ypost,2),size(ypost,3));
%         % Trial-averaging for AEF and other datasets
%         ypre_avg = mean(ypre,3);
%         ypost_avg = mean(ypost,3);
%
%         ypre_con = reshape(ypre,nc,size(ypre,2)*size(ypost,3));
%         ypost_con = reshape(ypost,nc,size(ypost,2)*size(ypost,3));
%
%
%         ypost = ypost_avg;

% Noise covariance initialization via VBFA
plotflag=1;
nl=5;
nem=100;
NM = 100;
%
% %         for AEF
% %         [b,lam,sig,yc,cy,bet,weight,mlike,xubar]=nut_reg_vbfa(ypre_avg,nl,nem,plotflag);
%       for picture-naming
[b,lam,sig,yc,cy,bet,weight,mlike,xubar]=nut_reg_vbfa(ypre,nl,nem,plotflag);


%         sigu0 = rand(size(sig));
%         sigu0 = 0.5 .* sig; % scale to 10-50%
sigu0 = sig; % scale to 10-50%

%         % Noise Learning Champagne with Convex Bounding
%         nem=100;
%         coup = 0;
%         noup = 3;
%         plot_on=1;
%         ncf = 0;  %% hom
%         vcs = 0;
%         sigu0 = rand(nc,nc);
% %         % for AEF
% %         [gamma,x,w,c,like,sigu_nlchamp]=champ_noise_up(double(ypre_avg),double(LF2),double(sigu0),nem,nd,vcs,plot_on,coup,noup,ncf);
%         % for picture-naming
%         [gamma,x,w,c,like,sigu_nlchamp]=champ_noise_up(double(ypre),double(LF2),double(sigu0),nem,nd,vcs,plot_on,coup,noup,ncf);
% %         sigu0 = 0.5 .* sigu_nlchamp;
%         sigu_add = rand(nc,nc);
%         sigu0 = sigu_nlchamp+sigu_add;

%% ---------------------------------------------
% Normalizing the leadfiled
% DW = sqrt(sum(L.^2));
% L = L*diag(1./DW);

%%
noise_cov = double(sigu0);

% noise_cov = rand(size(sigu0));
Y_total = double(ypost_reshape);
L = double(LF2);
NT = size(Y_total,2);
SigmaY = Y_total * Y_total'/NT;
DW = sqrt(sum(L.^2));
% L = L*diag(1./DW);


% New version of the code 
tic
[X_total_toeplitz]= ...
    Cov_estimator_NeuRIPS_realdata_ver2(Y_total,L,'noise_cov',noise_cov,'noise_update_mode','Diagonal','temporal_update_mode','Toeplitz','max_iters',50,'epsilon',1e-8,'print',1,'print_figures',1);
time_toeplitz = toc;

% % New version of the code 
% tic
% [X_total_toeplitz,X_trials_toeplitz,Sigma_toeplitz, Gamma_toeplitz, B_toeplitz,Gamma_error_toeplitz, B_error_toeplitz,f_loglikelihood_toeplitz]= ...
%     Cov_estimator_NeuRIPS_realdata_ver1(Y_total,L,noise_cov,'Adaptive-kernel', 'Toeplitz');
% time_toeplitz = toc;

% tic
% [X_total_toeplitz,X_trials_toeplitz,Sigma_toeplitz, Gamma_toeplitz, B_toeplitz,Gamma_error_toeplitz, B_error_toeplitz,f_loglikelihood_toeplitz]= ...
%     CovarianceEstimator_sensor_SigmaY_X_stable_realdata(Y_total,L,noise_cov,'sum_of_rank_one', 'Toeplitz');
% time_toeplitz = toc;

%         tic
%         [X_total_toeplitz,X_trials_toeplitz,Sigma_toeplitz, Gamma_toeplitz, B_toeplitz,Gamma_error_toeplitz, B_error_toeplitz,f_loglikelihood_toeplitz]= ...
%             CovarianceEstimator_sensor_SigmaY_X_stable(Y_total,L,SigmaY,Lambda_0,'Geodesic', 'Toeplitz');
%         time_toeplitz = toc;


X_total_toeplitz = diag(1./DW) * real(X_total_toeplitz);

test = input('stop');

%%
% tic
% [X_total_geodesic,X_trials_geodesic,Sigma_geodesic, Gamma_geodesic, B_geodesic,Gamma_error_geodesic, B_error_geodesic,f_loglikelihood_geodesic]= ...
%     CovarianceEstimator_sensor_SigmaY_X_stable_realdata(Y_total,L,SigmaY,noise_cov,'sum_of_rank_one', 'Geodesic');
% time_geodesic = toc;

%         tic
%         [X_total_geodesic,X_trials_geodesic,Sigma_geodesic, Gamma_geodesic, B_geodesic,Gamma_error_geodesic, B_error_geodesic,f_loglikelihood_geodesic]= ...
%             CovarianceEstimator_sensor_SigmaY_X_stable(Y_total,L,SigmaY,Lambda_0,'Geodesic', 'Geodesic');
%         time_geodesic = toc;

% X_total_geodesic = diag(1./DW) * real(X_total_geodesic);
% 
% corr_B_geodesic = corr(real(B_geodesic(:)),B_0(:));
% NMSE_B_geodesic = norm(B_0-real(B_geodesic),'fro')^2/norm(B_0,'fro')^2;
% 
% corr_Gamma_geodesic = corr(real(Gamma_geodesic(:)),Gamma_0(:));
% NMSE_Gamma_geodesic = norm(Gamma_0-real(Gamma_geodesic),'fro')^2/norm(Gamma_0,'fro')^2;
% 

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



%     end
% end


