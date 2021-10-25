clc
clear
close all

% Setting directoris.

folder = ['/data/research_meg3/yijing_gao/For_Ali/Dugh/NeurIPS-2021-spatiotemporal-Bayesian/Dugh_Realdata'];


%% Load Data

%nuts = load('AEF.mat'); % raw AEF data
nuts = load('VEF_1_70Hz.mat'); % filtered VEF data



rt = nuts.meg.srate;
data = nuts.meg.data;

[sp nc nt] = size(data);
for i=1:nc
    for j =1:nt
        data(:,i,j) = detrend(data(:,i,j));
    end
end

%%  Filter the data if it's raw
% data = filter_cc(data,'firls','bp','auto',1,70,rt,1);

%% set time point and preprocess leadfields
nsp = size(nuts.meg.latency,1);
pre_start=0.1*rt;
pre_end = find(0==nuts.meg.latency);

% for AEF task
% post_start = pre_end+.05*rt;
% post_end = pre_end+.15*rt;

% for VEF  task
post_start = pre_end+.025*1200;
post_end = pre_end+.275*1200;

[nc,nd,nv] = size(nuts.Lp);

LF = reshape(nuts.Lp,nc,nd*nv);
for i=1:size(LF,2)
    LF(:,i) = LF(:,i)./sqrt(sum(LF(:,i).^2));
end

LF_eLORETA_tmp = reshape(LF,nc,nd,nv);
LF_eLORETA = permute(LF_eLORETA_tmp,[1 3 2]);
LF2 = LF;

voxels = nuts.voxels;
voxelsize=nuts.voxelsize;
coreg=nuts.coreg;
bands=[1 70];
srate=nuts.meg.srate;
timepts = nuts.meg.latency(post_start:post_end,1);


%% Preprocessing & estimate noise variance
num_trials = [5 10 20 40 60 nt];
for i_num = 1:3
    nameresult = [folder, '/numoftrials_',  num2str(num_trials(i_num))];
    
    rand_sort = randperm(size(data,3));
     
    mkdir(nameresult);
    
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
    y_prestim = mean(ypre,3);
    y_poststim = ypost_reshape;

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
%     % for AEF
%     [b,lam,sig,yc,cy,bet,weight,mlike,xubar]=nut_reg_vbfa(ypre_avg,nl,nem,plotflag);
    % for VEF
    [b,lam,sig,yc,cy,bet,weight,mlike,xubar]=nut_reg_vbfa(y_prestim,nl,nem,plotflag);
    sigu0 = sig; 

    %% Preprocess inputs
    noise_cov = double(sigu0);
    Y_total = double(y_poststim);
    L = double(LF2);
    NT = size(Y_total,2);
    SigmaY = Y_total * Y_total'/NT;
    DW = sqrt(sum(L.^2));
        
  
    %% Dugh-T update mode
    tic
    [X_total_toeplitz]= ...
        Cov_estimator_NeuRIPS_investigate(Y_total,L,'noise_cov',noise_cov,'noise_update_mode','Diagonal','temporal_update_mode','Toeplitz','max_iters',50,'epsilon',1e-8,'print',1,'print_figures',1);
    time_toeplitz = toc;
    
    X_total_toeplitz = diag(1./DW) * real(X_total_toeplitz);
    
    % average trials for reconstructed sources (non-phase-locked data)
    X_total_toeplitz_reshape = reshape(X_total_toeplitz,size(X_total_toeplitz,1),post_end-post_start+1,num_trials(i_num));
    X_total_toeplitz_avg = mean(X_total_toeplitz_reshape,3);
    
    % save the results for Dugh-T
    x = X_total_toeplitz_avg;
    x = real(x);
    cd(nameresult);
    sa = [];
    sa{1}(:,:,1,1) = x(1:nd:end,:);
    sa{1}(:,:,1,2) = x(2:nd:end,:);
    sa{1}(:,:,1,3) = x(3:nd:end,:);
    
    save(['NeurIPS_Toeplitz_',num2str(num_trials(i_num)),'trials.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands','rand_sort')
    
    
    %% Dugh-G update mode
    
    % New version of the code
    tic
    [X_total_geodesic]= ...
        Cov_estimator_NeuRIPS_investigate(Y_total,L,'noise_cov',noise_cov,'noise_update_mode','Diagonal','temporal_update_mode','Geodesic','max_iters',50,'epsilon',1e-8,'print',1,'print_figures',1);
    time_geodesic = toc;
    
    X_total_geodesic = diag(1./DW) * real(X_total_geodesic);
    
    % average trials for reconstructed sources (non-phase-locked data)
    X_total_geodesic_reshape = reshape(X_total_geodesic,size(X_total_geodesic,1),post_end-post_start+1,num_trials(i_num));
    X_total_geodesic_avg = mean(X_total_geodesic_reshape,3);
    
    
    % save the results for Dugh-G
    x = X_total_geodesic_avg;
    x = real(x);
    cd(nameresult);
    sa = [];
    sa{1}(:,:,1,1) = x(1:nd:end,:);
    sa{1}(:,:,1,2) = x(2:nd:end,:);
    sa{1}(:,:,1,3) = x(3:nd:end,:);
    
    save(['NeurIPS_Geodesic_',num2str(num_trials(i_num)),'trials.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands','rand_sort')

    
    %% =============================================== %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ====Heteroscedastic Champagne from Chang==== %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Champagne with CC noise learning and Convex bounding
    ypost_champ = zeros(size(ypost,1),size(ypost,2),size(ypost,3));
    for g = 1:size(ypost,3)
        ypost_champ(:,:,g)=reshape(ypost_reshape(:,g),size(ypost,1),size(ypost,2));
    end

    ypost_champ_perm = permute(ypost_champ,[2 1 3]);
    Y_champ = reshape(ypost_champ_perm,size(ypost_champ_perm,1),[]);
    coup = 0;
    noup = 3;
    ncf = 1;  %% hetero
    NM = 100;
    tic
    [gamma,x_ChangHeteroChamp_ori,w,c,LK_ChampCC1,noise_hetero]=champ_noise_up(double(Y_champ),L,noise_cov,NM,nd,1,1,coup,noup,ncf);
    time_ChangHeteroChamp = toc;
    
    X_ChangHeteroChamp = diag(1./DW) * real(x_ChangHeteroChamp_ori);
    
    % average trials for reconstructed sources (non-phase-locked data)
    X_ChangHeteroChamp_reshape = reshape(X_ChangHeteroChamp,size(X_ChangHeteroChamp,1),post_end-post_start+1,num_trials(i_num));
    X_ChangHeteroChamp_avg = mean(X_ChangHeteroChamp_reshape,3);
    
    % save the results for ChangHeteroChamp
    x = X_ChangHeteroChamp_avg;
    x = real(x);
    cd(nameresult);
    sa = [];
    sa{1}(:,:,1,1) = x(1:nd:end,:);
    sa{1}(:,:,1,2) = x(2:nd:end,:);
    sa{1}(:,:,1,3) = x(3:nd:end,:);
    
    save(['NeurIPS_ChangHeteroChamp_',num2str(num_trials(i_num)),'trials.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands','rand_sort')
   
    %% =============================================== %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ========= eLORETA =========== %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    regu = 0.2;
    A = mkfilt_eloreta_v2(double(LF_eLORETA),regu);
    A = permute(A,[1 3 2]);
    A = reshape(A,size(A,1),size(A,2)*size(A,3));
    
    ypost_eLORETA = zeros(size(ypost,1),size(ypost,2),size(ypost,3));
    for g = 1:size(ypost,3)
        ypost_eLORETA(:,:,g)=reshape(ypost_reshape(:,g),size(ypost,1),size(ypost,2));
    end
    ypost_eLORETA_perm = permute(ypost_eLORETA,[2 1 3]);
    Y_eLORETA = reshape(ypost_eLORETA_perm,size(ypost_eLORETA_perm,1),[]);
    
    x_eLORETA = A' * Y_eLORETA;
    
    % average trials for reconstructed sources (non-phase-locked data)
    x_eLORETA_reshape = reshape(x_eLORETA,size(x_eLORETA,1),post_end-post_start+1,num_trials(i_num));
    x_eLORETA_avg = mean(x_eLORETA_reshape,3);
    
    % save the data
    x = x_eLORETA_avg;
    x = real(x);
    cd(nameresult);
    sa = [];
    sa{1}(:,:,1,1) = x(1:nd:end,:);
    sa{1}(:,:,1,2) = x(2:nd:end,:);
    sa{1}(:,:,1,3) = x(3:nd:end,:);
    
    save(['NeurIPS_eLORETA_',num2str(num_trials(i_num)),'trials.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands','rand_sort')
    
    
    
    %% =============================================== %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ========= sLORETA =========== %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LsLORETA = nuts.Lp;
    
    ypost_sLORETA = zeros(size(ypost,1),size(ypost,2),size(ypost,3));
    for g = 1:size(ypost,3)
        ypost_sLORETA(:,:,g)=reshape(ypost_reshape(:,g),size(ypost,1),size(ypost,2));
    end
    ypost_sLORETA_perm = permute(ypost_sLORETA,[2 1 3]);
    Y_sLORETA = reshape(ypost_sLORETA_perm,size(ypost_sLORETA_perm,1),[]);
    
    datasLORETA.Ryy= Y_sLORETA*Y_sLORETA'./size(Y_sLORETA,2);
    %     [weightsLORETA]= sLORETA(double(LsLORETA),data,'auto');
    [weightsLORETA]= sLORETA(double(LsLORETA),datasLORETA);
    wsLORETA = reshape(weightsLORETA,size(weightsLORETA,1),size(weightsLORETA,2)*size(weightsLORETA,3));
    x_sLORETA = wsLORETA'*Y_sLORETA;
    
    % average trials for reconstructed sources (non-phase-locked data)
    x_sLORETA_reshape = reshape(x_sLORETA,size(x_sLORETA,1),post_end-post_start+1,num_trials(i_num));
    x_sLORETA_avg = mean(x_sLORETA_reshape,3);
    
    % save the data
    x = x_sLORETA_avg;
    x = real(x);
    cd(nameresult);
    sa = [];
    sa{1}(:,:,1,1) = x(1:nd:end,:);
    sa{1}(:,:,1,2) = x(2:nd:end,:);
    sa{1}(:,:,1,3) = x(3:nd:end,:);
    
    save(['NeurIPS_sLORETA_',num2str(num_trials(i_num)),'trials.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands','rand_sort')
    
     %% MCE test
    ypost_MCE = zeros(size(ypost,1),size(ypost,2),size(ypost,3));
    for g = 1:size(ypost,3)
        ypost_MCE(:,:,g)=reshape(ypost_reshape(:,g),size(ypost,1),size(ypost,2));
    end
    
    ypost_MCE_perm = permute(ypost_MCE,[2 1 3]);
    Y_MCE = reshape(ypost_MCE_perm,size(ypost_MCE_perm,1),[]);
    
    plot_on=1;
   
    [x_MCE,w]=mce_ccai(Y_MCE,L,noise_cov,nem,nd,plot_on);
    x_MCE = real(x_MCE);
    
    % average trials for reconstructed sources (non-phase-locked data)
    x_MCE_reshape = reshape(x_MCE,size(x_MCE,1),post_end-post_start+1,num_trials(i_num));
    x_MCE_avg = mean(x_MCE_reshape,3);

    x = x_MCE_avg;
    cd(nameresult);
    sa = [];
    sa{1}(:,:,1,1) = x(1:nd:end,:);
    sa{1}(:,:,1,2) = x(2:nd:end,:);
    sa{1}(:,:,1,3) = x(3:nd:end,:);
    
    save(['NeurIPS_MxNE_',num2str(num_trials(i_num)),'trials.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands','rand_sort')

end

