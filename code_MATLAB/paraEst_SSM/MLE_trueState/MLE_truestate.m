function [mu1,cov1,sigma2] = MLE_truestate(U,ssmPar,semiEM)
% MLE from true state data
% Basic assumption: observe x_t at every time step
% %   x_{t+1} = x_t + dt*g(\theta,x_t)+ v_t,
%     Bn = x_{t+1} - x_t = Gn*\theta + v_t,


[d,tN] = size(U); 

stdF  = ssmPar.stdF; 
dt    = ssmPar.dt; 
K     = length(ssmPar.thetatrue);

% % % compute the covariance of theta likelihood
c1inv = zeros(K, K);
mu1   = zeros(K,1); 

Bt = zeros(d,tN-1); Gt = zeros(K,tN-1); 
for tt =1:tN-1
   [Bn,Gn] = BnGn(dt,U(:,tt),U(:,tt+1),K); 
   c1inv   = c1inv + Gn*Gn'; 
   mu1     = mu1 + Gn*Bn; 
   
   Bt(:,tt)= Bn;
   Gt(:,tt)= Gn; 
end
   c1inv = c1inv/(dt*stdF^2*tN);   % matrix A
   mu1   = mu1/(dt*stdF^2*tN);     % vector b

%% deal with singular A: regularize  --- pinv(A,tol) does this, and the default tol is better than tried numbers 
% mu0 =  pinv(c1inv)*mu1;  % the zero singular value will mu0 INFTY
%{
tol =10^(-3);  % tolerance in decompostion of Hessian singular values
[U,S,~] = svd(c1inv);     S    = sqrt(diag(S)); 
fprintf('\n Eigen-values of the likelihood matrix: \n');
disp(S');
if min(S)>tol
    mu1 = pinv(c1inv) *mu1;
    cov1= inv(c1inv);
else
    %     for k=3:5
    %         factor = (dt*stdF^2*tN)*10^(-k)
    %         fprintf('\n A regularization a*I is added with a= %3.6f \n', factor);
    %         cov0 = ( c1inv +factor*eye(size(c1inv)));
    %         mu0 = cov0\mu1; disp(mu0')
    %     end
    factor =0.5* S(end-1);
    fprintf('\n A regularization a*I is added with a= %3.7f \n', factor);
    cov1 = ( c1inv +factor*eye(size(c1inv)));
    mu1  = cov1\mu1; 
    % disp(mu1')
end
%}
    % regularMat =1/tN*eye(length(mu1)); 
   mu1   = pinv(c1inv)*mu1; 
   cov1  = c1inv; 
   
   vt    = Bt-mu1'*Gt; 
   sigma2= var(vt); 
   
if nargin ==3 && semiEM==1
    [mu1, ~] = change_semiEM(mu1,sigma2,dt,stdF); 
end

return


function [theta, theta2] = change_semiEM(mu1,sigma2,dt,stdF)  % TO BE DONE
% change of variables for the parameters of the semiEM alogrithm
c      = sqrt(dt)*stdF/sqrt(sigma2); 
theta2 = (1-c)/dt;  
theta = mu1*c;

diff = theta(2) -theta2; 
% fprintf('SemiEM: Difference theta2Euler -theta2semiEM: %2.4f\n', diff);

return