

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
%
% #s:
% nk = # of sensors
% nv = # of voxels
% nt = # of data time points


function [gamma,x,w,c,Gamma_error,l]=awsm_champ(y,f,sigu,nem,nd,vcs,plot_on);

if vcs==2 && nd>1
    [gamma x w l v c]=champ_mat(y,f,sigu,nem,nd,plot_on);
else
    [gamma x w l v c Gamma_error]=champ_vec(y,f,sigu,nem,nd,vcs,plot_on);
end

return




function [gamma,x,w,like,vvec,c,Gamma_error]=champ_vec(y,f,sigu,nem,nd,vcs,plot_on);

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

Gamma_old = diag(vvec); 
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
    
%    like(iem)=-.5*(sum(log(d))+nk*log(2*pi))-.5*sum(sum(y.*(invc*y)))/nt; 
%     like(iem)=-.5*(sum(log(max(d,eps1)))+nk*log(2*pi))-.5*sum(sum(invc.*cyy));  %4.77 4.18 y's distribution directly
    like(iem)= sum(log(max(d,eps1)))+sum(sum(invc.*cyy));  %4.77 4.18 y's distribution directly
    
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

    v=sum(reshape(vvec,nd,nv),1);   %voxel power
    
        
    if(plot_on)
    subplot(2,2,2);plot((1:nv),v);
    title(['Voxel power: ' num2str(nv) ' / ' num2str(nv)]);
    xlabel('voxel index');
    set(gca(),'XLim',[1 nv]);
    drawnow
    end
    
    Gamma = diag(vvec); 
    Gamma_error(iem) = norm(Gamma-Gamma_old,'fro'); 
    Gamma_old = Gamma; 
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




function [gamma,x,w,like,vvec,c]=champ_mat(y,f,sigu,nem,nd,plot_on,nupd);

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
    w=vmat*fc;
%    x=w*y;
%    x2=mean(x.^2,2);
%    z=sum(fc.*f',2);

    for iv=1:nv
        jv=(iv-1)*nd+1:iv*nd;             
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