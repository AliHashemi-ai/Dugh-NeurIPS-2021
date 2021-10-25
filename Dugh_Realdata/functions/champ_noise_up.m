

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learn voxel covariances in Champagne
% 
% � 2011 Convex Imaging
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Output: 
% gamma(nd,nd,nv) = voxel covariance matrices
% x(nv*nd,nt) = voxel current density 
% w(nv*nd,nk) = reconstruction filter
%
% Input:
% y(nk,nt) = sensor data 
% f(nk,nv*nd) = lead field matrix
% sigu(nk,nk) = noise covariance matrix
% nd = # of orientations
% vcs = voxel covariance structure: 0 = scalar, 1 = diagonal, 2 = general 
% nem = maximum # of iterations
% noup = 0, no noise update, 1 ML, 2 MK, 3 Convex, 4 EM
%% add two more parameter for input
% coup: covaraince update methods, 0=covex boounding, 1 = EM, 2 = Mackay
% noup: noise update methods, 0 = no update, 1 = EM, 2 = Ali.  3=convex,
% 4=em
% ncf: noise covavariance form, 0 = scalar, 1 = heter, 2 = full
%
% #s:
% nk = # of sensors
% nv = # of voxels
% nt = # of data time points




function [gamma,x,w,c,l, sigu]=champ_noise_up(y,f,sigu,nem,nd,vcs,plot_on,coup,noup,ncf);

if vcs==2 && nd>1
    [gamma x w l v c sigu]=champ_mat(y,f,sigu,nem,nd,plot_on,coup,noup,ncf);
else
    [gamma x w l v c sigu]=champ_vec(y,f,sigu,nem,nd,vcs,plot_on,coup,noup,ncf);
end

return




function [gamma,x,w,like,vvec,c,sigu]=champ_vec(y,f,sigu,nem,nd,vcs,plot_on,coup,noup,ncf);

eps1=1e-8;
[nk nvd]=size(f);
nv=nvd/nd;
nt=size(y,2);

cyy=y*y'/nt;

% Initialize voxel variances

f2=sum(f.^2,1);             %which filter used to estimated 
invf2=zeros(1,nvd);
ff=find(f2>0);
invf2(ff)=1./f2(ff);
f=double(f);
w=spdiags(invf2',0,nvd,nvd)*f';
%w=spdiags(1./sum(f.^2,1)',0,nvd,nvd)*f';
inu0=mean(mean((w*y).^2));
%gamma=inu0*repmat(eye(nd,nd),[1,1,nv]);
vvec=inu0*ones(nvd,1);

insigu = inv(sigu);
if sigu==0
    insigu = 0;
end
% Learn voxel variances
if(plot_on)
figure;
end

like=zeros(nem,1);
for iem=1:nem
%    iem
    vmat=spdiags(vvec,0,nvd,nvd);
    c=f*vmat*f'+sigu;   %used for 4.74
%    [p d q]=svd(c);
    [p d]=eig(c);%used for 4.74
    d=max(real(diag(d)),0);%used for 4.74
    invd=zeros(nk,1);%used for 4.74
    ff=find(d>=eps1);%used for 4.74
    invd(ff)=1./d(ff);   %concave function parameters 1/a1,1/a2,....1/aN %used for 4.74
    invc=p*spdiags(invd,0,nk,nk)*p';        %model data convariance y %used for 4.74
    invc = inv(c);
    
%     [p_s d_s]=eig(sigu);%used for 4.74
%     d_s=max(real(diag(d_s)),0);%used for 4.74
%     invd_s=zeros(nk,1);%used for 4.74
%     ff=find(d_s>=eps1);%used for 4.74
%     invd_s(ff)=1./d_s(ff);   %concave function parameters 1/a1,1/a2,....1/aN %used for 4.74
%     invsigu=p_s*spdiags(invd_s,0,nk,nk)*p_s';        %model data convariance y %used for 4.74
%     %invc = inv(c);
%     
%    like(iem)=-.5*(sum(log(d))+nk*log(2*pi))-.5*sum(sum(y.*(invc*y)))/nt; 
    like(iem)=-.5*(sum(log(max(d,eps1)))+nk*log(2*pi))-.5*sum(sum(invc.*cyy));  %4.77 4.18 y's distribution directly
    
    
    
    
    
    if(plot_on)
    subplot(2,2,1);plot((1:iem),like(1:iem));
    title(['Likelihood: ' int2str(iem) ' / ' int2str(nem)]);
    xlabel('iteration');
    set(gca(),'XLim',[0 iem]);
    end
    
    fc=f'*invc;         %4.82
    w=vmat*fc;          %4.82 update s
    x=w*y;              %4.82 source distribution
    x2=mean(x.^2,2);    %xk distribution is equal to the posterior mean of the source distribution
    z=sum(fc.*f',2);    %update Z 4.79
    
    if 0 == coup  %% convex update
        if vcs==0
            x20=sum(reshape(x2,nd,nv),1);
            z0=sum(reshape(z,nd,nv),1);
            v0=zeros(size(z0));
            ff=find(z0>0);
            v0(ff)=sqrt(max(x20(ff)./z0(ff),0)); % CN 10/2012 added max,0
    %        v0(ff)=sqrt(x20(ff)./z0(ff));
            vvec=reshape(ones(nd,1)*v0,nvd,1);
        else
            vvec=zeros(size(x2));
            ff=find(z>0);
            vvec(ff)=sqrt(max(x2(ff)./z(ff),0)); % CN 10/2012 power 4.88 update v
    %                vvec(ff)=sqrt(x2(ff)./z(ff));
        end

        vvec = real(vvec);
        v=sum(reshape(vvec,nd,nv),1);   %voxel power

        
    elseif 1 == coup %% EM update
        
        igam = sum((eye(nvd)-w*f).*vmat,2);
        vvec = (igam)+x2;
        v = sum(reshape(vvec,nd,nv),1);   %voxel power

        
    elseif 2 == coup %% Mackay
         
       ZZ = sum(f'*invc.*f',2);%insigu*invcovx*f'*f;

        dzz = vvec.*ZZ;
        ff = find(dzz>0);
        vvec(ff) = x2(ff)./dzz(ff);
        v = sum(reshape(vvec,nd,nv),1);   %voxel power

    end
    
    
    
%% noise update noup=1 is EM; noup=2 is Ali; noup=3 is Mackay
   
    
    if  1==noup  %Ali's update rules for noise 
       	n0 =y-f*x;
        mn = mean(n0.^2,2);
        
        ff = find(vvec>0);
        invvv = vvec; %zeros(size(vvec));
        invvv(ff) = 1./vvec(ff);
        
        covx = f'*insigu*f+diag(invvv);%inv(diag(vvec));

        invcovx = inv(covx);
        Reg = nk-nvd+sum(diag(invcovx).*invvv);
        
        sigu = (sum(mn)./Reg);
        
        sigu = eye(nk)*sigu;
        
    elseif 2==noup  %% Mackay
        if ncf ==0  % scalar
            fw=eye(nk)-f*w;
            sigy1=mean(sum((fw*cyy).*fw,2));  
            fvfinsig = mean(sum(f*vmat*f'.*inv(sigu),2))+1;
            sigu =fvfinsig*sigy1*eye(nk);
        elseif ncf == 1 % vector
            fw=eye(nk)-f*w;
            sigy1=sum((fw*cyy).*fw,2);
            % fvfinsig = sum(f*vmat*f'.*invsigu,2)+1;
            fvfinsig = sum(f*vmat*f'.*inv(sigu),2)+1;
            sigu000 = sigu;
            sigu = diag(fvfinsig.*sigy1); 
            norm(sigu-sigu000);
            
            
%             invgam = inv(f'*inv(sigu)*f+inv(vmat));
%             term1 = sum(pinv(invgam'*f')*vmat.*f,2);
%             sigu =diag(term1.*sigy1);
        end
        
    elseif 3 == noup %% convex based
        if ncf ==0  % scalar
            fw=eye(nk)-f*w;
            sigy1=mean(sum((fw*cyy).*fw,2));   
            cc = 1./mean(diag(invc));
            sigu =sqrt(cc*sigy1)*eye(nk);
        elseif ncf == 1 % vector
            fw=eye(nk)-f*w;
            sigy1=(sum((fw*cyy).*fw,2));
            cc = 1./diag(invc);
            sigu000 = sigu;
            sigu = diag(sqrt(cc.*sigy1));
            norm(sigu-sigu000);
            
           
         end
     
    elseif 4 == noup %% bayesian update rule EM
        if ncf == 0 %% scalar
            fw=eye(nk)-f*w;
            sigy1=mean(sum((fw*cyy).*fw,2));

            fgf=sum((fw*f*vmat).*f,2);
            ilam=sigy1+mean(fgf);  %% ilam = sigy1+fgf; for diagnol
            sigu=ilam*eye(nk);
        elseif ncf == 1 %% heter
            y1=y-f*x;  
            sigy1 = diag(sum(y1.*y1,2)/nt);

            fw=eye(nk)-f*w;

            fgf=diag(sum((fw*f*vmat).*f,2));

            sigu0 = sigy1+(fgf);
            sigu = double(sigu0);
        elseif ncf == 2 %% matrix
            
            y1=y-f*x;  
            sigy1 = (y1*y1')/nt;

            fw=eye(nk)-f*w;

            fgf=fw*f*vmat*f';%diag(sum((fw*f*vmat).*f,2));

            sigu0 = sigy1+(fgf);
            sigu = double(sigu0);
            
            %sigu = sigu./sum(diag(sigu));
            
        end
%     elseif 5 == noup
%         if ncf ==0  % scalar
%             fw=eye(nk)-f*w;
%             sigy1=mean(sum((fw*cyy).*fw,2));   
%             cc = 1./mean(diag(invc));
%             sigu =sqrt(cc*sigy1)*eye(nk);
%         elseif ncf == 1 % vector
%             fw=eye(nk)-f*w;
%             sigy1=(sum((fw*cyy).*fw,2));
%             cc = 1./diag(invc);
%             sigu = diag(sqrt(cc.*sigy1));
%         elseif ncf ==2 % full matrix
%             y1=y-f*x;  
%             sigy1 = y1*y1'/nt;
%             cc = c*sigy1;
%             [p d q] = svd(cc);
%             d =diag(d);
%             d(21:end) = 0;
%             d = sqrt(d);
%             sigu = (p*diag(d)*q');
%             %sigu = (sqrtm(c'*sigy1));
%         end
    end
        
    
    if(plot_on)
    subplot(2,2,2);plot((1:nv),v);
    title(['Voxel power: ' num2str(nv) ' / ' num2str(nv)]);
    xlabel('voxel index');
    set(gca(),'XLim',[1 nv]);
    drawnow
    end
    
%    lam=inv(sigu);nu=inv(vmat);
%    gam=f'*lam*f+nu;
%    igam=inv(gam);
%    w1=igam*f'*lam;
%    x1=w1*y;
%    x2=mean(x1.^2,2)+diag(igam);
%    disp([max(max(abs(w1-w))) max(abs(x2-v))]);   
end
   
if nd==1
    gamma=reshape(vvec,1,1,nv);
else
    gamma=zeros(nd,nd,nv);
    for iv=1:nv
        gamma(:,:,iv)=diag(vvec((iv-1)*nd+1:iv*nd));
    end
end

return




function [gamma,x,w,like,vvec,c,sigu]=champ_mat(y,f,sigu,nem,nd,plot_on,coup,noup);

eps1=1e-8;
[nk nvd]=size(f);
nv=nvd/nd;
nt=size(y,2);

cyy=y*y'/nt;

% Initialize voxel variances

f2=sum(f.^2,1);
invf2=zeros(1,nvd);
ff=find(f2>0);
invf2(ff)=1./f2(ff);
w=spdiags(invf2',0,nvd,nvd)*f';
%w=spdiags(1./sum(f.^2,1)',0,nvd,nvd)*f';
inu0=mean(mean((w*y).^2));
%gamma=inu0*repmat(eye(nd,nd),[1,1,nv]);
v=zeros(nv,1);
vmat=double(inu0)*speye(nvd,nvd);

% Learn voxel variances
 if(plot_on)
figure;
 end
 
like=zeros(nem,1);
for iem=1:nem
%     iem
%    disp(full(vmat));pause
%    vmat=spdiags(v,0,nvd,nvd);
    c=f*vmat*f'+sigu;
%    [p d q]=svd(c);
    [p d]=eig(c);
    d=max(real(diag(d)),0);
    invd=zeros(nk,1);
    ff=find(d>=eps1);
    invd(ff)=1./d(ff);
    invc=p*spdiags(invd,0,nk,nk)*p';

%    like(iem)=-.5*(sum(log(d))+nk*log(2*pi))-.5*sum(sum(y.*(invc*y)))/nt;
    like(iem)=-.5*(sum(log(max(d,eps1)))+nk*log(2*pi))-.5*sum(sum(invc.*cyy));  
    
    if(plot_on)
    subplot(2,2,1);plot((1:iem),like(1:iem));
    title(['Likelihood: ' int2str(iem) ' / ' int2str(nem)]);
    xlabel('iteration');
    set(gca(),'XLim',[0 iem]);
    end
    fc=f'*invc;
    w=real(vmat*fc);
    
%    x=w*y;
%    x2=mean(x.^2,2);
%    z=sum(fc.*f',2);

    for iv=1:nv
        jv=((iv-1)*nd+1:iv*nd);
%        x2=x(jv,:)*x(jv,:)'/nt;
        x2=w(jv,:)*cyy*w(jv,:)';
        z=fc(jv,:)*f(:,jv);  %4.58 update Z
        
        [pz dz]=eig(z);
        dz5=sqrt(max(real(diag(dz)),0));
%        dz5=sqrt(abs(diag(dz)));
        invdz5=zeros(nd,1);
        ff=find(dz5>=eps1);
        invdz5(ff)=1./dz5(ff);
        z5=pz*diag(dz5)*pz';
        invz5=pz*diag(invdz5)*pz';

        [px dx]=eig(z5*x2*z5);
        dx5=sqrt(max(real(diag(dx)),0));
%        dx5=sqrt(abs(diag(dx)));
        cx5=px*diag(dx5)*px';
        vmat(jv,jv)=invz5*cx5*invz5;
        v(iv)=sum(diag(vmat(jv,jv)));
    end
if(plot_on)
    subplot(2,2,2);plot((1:nv),v);
    title(['Voxel power: ' num2str(nv) ' / ' num2str(nv)]);
    xlabel('voxel index');
    set(gca(),'XLim',[1 nv]);
    drawnow
end
%    lam=inv(sigu);nu=inv(vmat);
%    gam=f'*lam*f+nu;
%    igam=inv(gam);
%    w1=igam*f'*lam;
%    x1=w1*y;
%    x2=mean(x1.^2,2)+diag(igam);
%    disp([max(max(abs(w1-w))) max(abs(x2-v))]);   
end

x=w*y;
gamma=zeros(nd,nd,nv);
for iv=1:nv
    jv=(iv-1)*nd+1:iv*nd;
    gamma(:,:,iv)=vmat(jv,jv);
end
vvec=diag(vmat);

return


  
%                sigu00 = y1*y1'./size(y1,2);
%                sigu0 = diag(mean(y1*y1',2))-mean(y1,2)*mean(y1,2)';
%                sigu000 = double(sigu0);
%                sigu0 = max(eig(sigu00'))*eye(size(sigu00));

%                sigu = double(max(eig(sig'))*eye(size(sig)));
               
%                sigy1=mean(mean(y1.^2));
%                
%         fw=eye(nk)-f*w;
%         sigy1=mean(sum((fw*cyy).*fw,2));
%         fgf=v*sum((fw*f).*f,2);
%         ilam=sigy1+mean(fgf);
%         sigu=ilam*eye(nk);