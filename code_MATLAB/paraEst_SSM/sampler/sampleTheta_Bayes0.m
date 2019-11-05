%-------------------------------------------------------------------
function [thetaP,nIter] = sampleTheta_Bayes(U,ssmPar,K,dt, prior)
% sample the posterior of theta; 
% Input
%    U    - a trajectory of the state
%    K    - dimension of parameter theta
%    dt   - time step size
%   prior - prior of theta
% 
% Because of linear dependence on theta, the likelihoood is computed
% analytically


if nargin <5; prior.flag = 0; end   % if no prior, use N(0,I)  

[~,tN] = size(U); 

% % % compute the covariance of theta likelihood
c1inv = zeros(K, K);      % covariance of likelihood. 
mu1   = zeros(K,1); 

stdF  = ssmPar.stdF; 

for tt =1:tN-1
    U0     = U(:,tt); U1 = U(:,tt+1);
   [Bn,Gn] = BnGn(dt,U0,U1,K); 
   c1inv   = c1inv + Gn*Gn'; 
   mu1     = mu1 + Gn*Bn; 
end
   c1inv = c1inv/(dt*stdF^2*tN); mu1= mu1/(dt*stdF^2*tN);

%%{
 %% regularize singular c1inv; ==== Not helping; pinv(A,tol) does this, and default tol is better 
 tol =1e-3;  % tolerance in decompostion of Hessian singular values
[U,S,~] = svd(c1inv);     S    = sqrt(diag(S));
if min(S)>tol
    mu1 = pinv(c1inv) *mu1;  Cov1_chol= chol(c1inv);
else              
    fprintf('\n Singular values of the likelihood matrix: \n'); 
    disp(S'); fprintf('Progress:    ');
    ii   = find(abs(S)>tol);   Sinv = 0*S + 1/tol;
    Sinv(ii)= 1./S(ii);           Cov1 = U'*diag(Sinv)*U;
    Cov1_chol= chol(Cov1);
    mu1 = Cov1*mu1 ; % c1inv\ mu1;   
end
%}

%% posterior = likelihood and prior 
% % % % combine likelihood with prior
if prior.flag ==0    % use Gaussian prior
   covPriorInv = diag(1./prior.sigma); 
   covInv = inv(covPriorInv+c1inv);  covInv_chol = chol(covInv);
   muPost = (covPriorInv +c1inv) \ (c1inv*mu1+ covPriorInv*prior.mu);
   theta  = randn(K,1);
   thetaP = muPost+ covInv_chol*theta;

elseif prior.flag==1  % use prior:[log(N(m0,s0^2)),N(m123,s123^2), -logN(m4,s4^2)]
   % use the prior as proposal + the rejection method 
   nIter = 0;  maxN= 10000;
   % compute the covariance  to be used in mvnpdf
   while nIter < maxN+1
       temp    = randn(K,1).*sqrt(prior.sigma) + prior.mu;
       thetaP  = [exp(temp(1)); temp(2:4); -exp(temp(5))];
       ratio   = mvnpdf(thetaP, mu1, Cov1); % *det(Cov1);
       if ratio < 10^(-3)
           fprintf('Posterior sampling failure: Lognormal\n');
           break;
       end
       if rand(1) <ratio; disp('accept');break; end
       nIter = nIter + 1;
       if nIter== maxN
           fprintf('Maximum iteration in A-R sampling. MLE used \n');
           thetaP = mu1; % use the MLE if failure of Acceptance-Rejction
           thetaP(1) = abs(mu1(1));  thetaP(5) = -abs(mu1(5));
           break;
       end
   end
elseif prior.flag == 2   % uniform prior: sample from truncated normal (bc. likelihood is normal)
     thetaP = IS_Unifprior(prior,mu1,c1inv);  % use importance sampling 
     %%%   use truncated Gaussian 
%    lb   = Cov1_chol\(prior.lb- mu1); ub =Cov1_chol\(prior.ub -mu1);
%    temp = trandn(lb,ub);
%    thetaP = mu1+ Cov1_chol*temp;
   if thetaP(1)<0 
       fprintf('theta0  wrong sign\n');
  %     keyboard;
   end
end

end



function sampleoutput = IS_Unifprior(prior,mu,cov)
% importance sampling for the case with Uniform prior
lb     = prior.lb; ub = prior.ub; 
K      = length(lb); 
numS   = 10; 
sample = zeros(K,numS); wts = ones(1,numS); 
for n = 1:numS 
    temp   = rand(K,1).*(ub-lb) + lb;
    sample(:,n) = temp; 
    sample2 = temp - mu;
    temp    = pinv(cov)*sample2;
    wts(n)  = sum(sample2.*temp);
end
    wts     = wts - max(wts);
    wts     = exp(wts);
    wts     = wts/sum(wts);

indres = myResampling(wts); 
sample = sample(:,indres); 
sampleoutput = sample(:,1); 
end


