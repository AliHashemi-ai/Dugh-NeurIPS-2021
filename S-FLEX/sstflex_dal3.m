function [S, out] = sstflex_dal3(X, LL, para);
%
% Sparse Spatio-Temporal Basis Field Expansion (SST-FLEX), 
% Haufe et al. (2008)
%
% Synopsis:
%   S = sstflex(X, LL, R, para)
%   
% Arguments:
%   X : [M K] measurement matrix
%   LL: [M N NDUM] leadfield tensor
%   para
%     .eps: scalar between 0 (default, perfect fit) and 1 specifying the
%     desired fit
%     .B : [N K] matrix of basis function evaluations at positions .R
%     .R : [N NDUM] matrix of dipole positions. Needed if no basis .B is 
%     given. Then Gaussian basis functions centered at .R are taken. 
%
%   
% Returns:
%   S: [N NDUM K] estimated dipole moments
%   C: [NB NDUM*K] sparse estimated coefficients
%
% Stefan Haufe, 2007, 2008
%
%
% License
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see http://www.gnu.org/licenses/.

tic

[M N NDUM] = size(LL);

if nargin < 3
  para=[];
end

maxtol = 1e-8;
minresnorm = 1e-5;
if ~isfield(para, 'eps')
    para.eps = minresnorm;
else
    if para.eps < minresnorm 
        para.eps = minresnorm;
    end
    if para.eps > 1
        para.eps = 1;
    end
end

if ~isfield(para, 'factor')
    para.factor = 3/4;
end

if ~isfield(para, 'cv')
    para.cv = 0;
end

if ~isfield(para, 'eta')
    para.eta = 100;
end

if ~isfield(para, 'solver')
    para.solver = 'qn';
end

if para.cv == 1
    para.cv = 0;
end

if ~isfield(para, 'regpath')
    para.regpath = 0;
end

if para.cv
    para.regpath = 0;
end

[M_, K] = size(X);

if ~isfield(para, 'car')
    para.car = 1;    
end

if ~isfield(para, 'maxgroups')
    para.maxgroups = inf;    
end

if ~isfield(para, 'noisecov')
    para.noisecov = eye(M_);
else
    para.noisecov = para.noisecov / norm(para.noisecov);
end

if ~isfield(para, 'depthweighting')
    para.depthweighting = 'sLORETA';    
end

if para.car == 1
    H = eye(M_)-ones(M_) ./ M_;
    [UH SH VH] = svd(H);
    out.H =  SH(1:end-1, :)*VH';
else
    out.H = eye(M_);
end

out.spatfilt = real((out.H*para.noisecov*out.H')^(-0.5))*out.H;
out.backfilt = out.H'*real((out.H*para.noisecov*out.H')^(0.5));
    
M = size(out.spatfilt, 1);

fprintf('S-FLEX: preparing...');

out.L = out.spatfilt*reshape(LL, M_, N*NDUM);

if isequal(para.depthweighting, 'sLORETA')
  W = sloreta_invweights(reshape(out.L, M, N, NDUM));
  W = W / norm(W, 'fro');
else
  gamma = 0.7;
  
  L_temp = reshape(out.L, M, N, NDUM);
  W = sqrt((norms(L_temp(:, :, 1), 2, 1).^2+norms(L_temp(:, :, 2), 2, 1).^2 ...
    + norms(L_temp(:, :, 3), 2, 1).^2).^gamma)';
  W = spdiags(reshape(repmat(W, 1, 3), [], 1), 0, N*NDUM, N*NDUM);
  W = W / norm(W, 'fro');
end
out.L = out.L*W;
out.X = out.spatfilt*X;

compl = 0;
if ~isreal(out.X)
    compl = 1;
    out.X = [real(out.X) imag(out.X)];
    K = 2*K;
end

if para.cv
    FM = ceil(M/para.cv);
    per = randperm(M);
    for icv = 1:para.cv
        M_out{icv} = per((FM*(icv-1)+1):min(FM*(icv), M));
        M_in{icv} = setdiff(1:M, M_out{icv});
    end
end
M_out{para.cv+1} = [];
M_in{para.cv+1} = 1:M;

if isfield(para, 'B')
  NB = size(para.B, 2);   
  A = [];
  for ii = 1:NDUM
    A = [A out.L(:, (ii-1)*N+[1:N])*para.B];
  end   
else
%   try
    if ~isfield(para, 'sigma')
      para.sigma = 0.5:0.25:1.5;
    end  
    NB = length(para.sigma)*N;
    A = [];
    for isig = 1:length(para.sigma)
      H = gaussrbf_sparse(para.R', para.sigma(isig));
      H2 = [];
      for ii = 1:NDUM
        H2 = cat(3, H2, out.L(:, (ii-1)*N+[1:N])*H + out.L(:, (ii-1)*N+[1:N])*H');
      end 
      A = cat(2, A, H2);
    end
    clear H2;
    if isfield(para, 'sel')
        A = A(:, para.sel, :);
        NB_ = NB;
        NB = length(para.sel);
    end
    A = reshape(permute(A, [1 3 2]), [M NB*NDUM]);
%   catch
%       error('Either para.B or para.R needed.')
%   end
end

clear H LL X;
inds = {};
for ii = 1:NB
    inds{ii} = ([1:NDUM]-1)*NB+ii;
end
inds2 = reshape(reshape(1:NDUM*NB, [], NDUM)', [], 1);

if ~isfield(para, 'startlambda')
    para.startlambda = (M/2)*sqrt(K)/norm(out.X, 'fro');
end

fprintf('done.');
min_lambda = 1e-4;
out.opt_lambda = [];
out.opt_ind = 1;
out.trainerr = [];
out.testerr = [];
for icv = 1:para.cv+1
    y = reshape(out.X(M_in{icv}, :), [], 1);
    lambda = para.startlambda;
    factor = para.factor;
    summand = 0;
    kk = 0;

    if para.cv > 1
        if icv <= para.cv
            fprintf('\nOptimizing fold %d/%d', icv, para.cv);
        else
            fprintf('\nOptimizing final model');
        end
    else
        fprintf('\nOptimizing');
    end
    
    goup = 0;
    init = zeros(NDUM*K, NB);
    while true        
        AA = {@xforth, @xback, length(M_in{icv})*K, NB*NDUM*K};
%         tic
        if para.cv > 1 & icv > para.cv
            lambda = lambda*(para.cv/(para.cv-1));
        end
        [init,status]=dalsqgl(init, AA, y, lambda, 'eta', para.eta, 'solver', para.solver);
        init = full(init);
%         toc
        approx = A(M_in{icv}, :)*reshape(permute(reshape(init, NDUM, K, NB), [1 3 2]), NDUM*NB, K);
        
        noin = find(norms(reshape(permute(reshape(init, NDUM, K, NB), [1 3 2]), NDUM*NB, K), 2, 2));
        if ~isempty(noin)
%             approx = A(M_in{icv}, noin)*pinv(A(M_in{icv}, noin))*out.X(M_in{icv}, :);
            res = out.X(M_in{icv}, :)-approx;
            resnorm = norm(res, 'fro')/norm(out.X(M_in{icv}, :), 'fro');
        else
            resnorm = 1;
        end
%         fprintf('GOF = %8.5e\n', resnorm);

        kk = kk + 1;
        
        if para.cv 
            if icv < para.cv + 1
                out.trainerr(icv, kk) = norm(res, 'fro')^2;
                out.testerr(icv, kk) = norm(out.X(M_out{icv}, :)-A(M_out{icv}, :)*reshape(permute(reshape(init, NDUM, K, NB), [1 3 2]), NDUM*NB, K), 'fro')^2;
            else
                C = reshape(init', [], 1);
                out.lcurve = [sum(sqrt(sum(reshape(init, [], NB).^2))); resnorm];
                if ~compl
                    out.approx = approx;
                else
                    out.approx = [approx(:, 1:end/2) + sqrt(-1)*approx(:, ((end/2)+1):end)];
                end
                break
            end
        else
            if para.regpath
                C(:, kk) = reshape(init', [], 1);
                out.lcurve(:, kk) = [sum(sqrt(sum(reshape(init, [], NB).^2))); resnorm];
                if ~compl
                    out.approx(:, :, kk) = approx;
                else
                    out.approx(:, :, kk) = [approx(:, 1:end/2) + sqrt(-1)*approx(:, ((end/2)+1):end)];
                end
            else
                C(:, mod(kk, 2)+1) = reshape(init', [], 1);
                out.lcurve(:, mod(kk, 2)+1) = [sum(sqrt(sum(init.^2))); resnorm];
                if ~compl
                    out.approx(:, :, mod(kk, 2)+1) = approx;
                else
                    out.approx(:, :, mod(kk, 2)+1) = [approx(:, 1:end/2) + sqrt(-1)*approx(:, ((end/2)+1):end)];
                end
            end
        end
        
        if kk == 1
            if (resnorm-para.eps <= maxtol)
                goup = 1;
                factor = 1/factor;
            end
        else
            if (goup & para.eps-resnorm <= maxtol) | (~goup & resnorm-para.eps <= maxtol) | (~goup & length(find(status.info.spec)) > para.maxgroups) 
                break
            end
        end
        lambda = lambda*factor;
    end
    fprintf('\n');
    
    if para.cv 
        if icv == para.cv
            out.testerr(out.testerr < eps) = nan;
            out.te = mean(out.testerr);
            out.se = std(out.testerr)./sqrt(para.cv);
            [mi, out.opt_ind] = min(out.te);
%             out.opt_lambda = para.startlambda/factor.^(out.opt_ind-1);
            para.startlambda = para.startlambda*factor.^(out.opt_ind-1);
        end
    else
        [mi out.opt_ind] = min(abs(out.lcurve(2, :)-para.eps));
    end
    if ~para.regpath
        out.approx = out.approx(:, :, out.opt_ind);
        out.lcurve = out.lcurve(:, out.opt_ind);
        C = C(:, out.opt_ind);
    end
end

clear A

kk = size(C, 2);
C = reshape(C, [], NDUM*K*kk);
out.sel = find(norms(C, 2, 2) > 1e-6);

if isfield(para, 'B')
	S = reshape(full(W*reshape(para.B*C, [], K*kk)), N, NDUM, K, kk);
else
    if isfield(para, 'sel')
       C_ = C;
       C = spalloc(NB_, NDUM*K*kk, nnz(C_));
       C(para.sel, :) = C_;
       clear C_;
    end
  S = zeros(N, NDUM*K*kk);
  for isig = 1:length(para.sigma)
    H = gaussrbf_sparse(para.R', para.sigma(isig));
    S = S + H*C((isig-1)*N+[1:N], :);
	S = S + H'*C((isig-1)*N+[1:N], :);
  end
  clear C
  S = reshape(W*reshape(S, NDUM*N, K*kk), N, NDUM, K, kk);
  if compl
      S = S(:, :, 1:end/2, :) + sqrt(-1)*S(:, :, ((end/2)+1):end, :);
      out.X = out.X(:, 1:end/2) + sqrt(-1)*out.X(:, ((end/2)+1):end);
  end    
end

out.para = para;

toc

function xfo = xforth(x)
  Q = size(x, 2);
  [in1 indum] = find(x);
  in2 = unique(ceil(in1./(K)));
  L = length(in2)/NDUM;
  if L == 0
      xfo = zeros(length(M_in{icv})*K, Q);
  else
      xfo = reshape((A(M_in{icv}, in2)*reshape(permute(reshape(full(x(in1, :)), NDUM, K, L, Q), [1 3 2 4]), NDUM*L, K*Q)), length(M_in{icv})*K, Q);
  end
end
  
function xba = xback(x)
  xba = reshape(permute(reshape(A(M_in{icv}, :)'*reshape(x, length(M_in{icv}), K), NDUM, NB, K), [1 3 2]), NDUM*K*NB, []);
end

end
