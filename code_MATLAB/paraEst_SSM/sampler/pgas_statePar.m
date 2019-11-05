function [X,ess,theta] = pgas_statePar(numMCMC,obs,t0t1,ssmPar,Np,dt,...
                         tN,theta0,varObs, prior,sp,semiEM,ProgressOn)
% Runs the PGAS algorithm to estimate states/parameters of the state space model: 
%        x_{t+1} = x_t + g(\theta,x_t)+ v_t,
%         y_t     = x_t + W_t,
% Input
% % numMCMC     - number steps in MCMC
% % obs         - observations:   size =[ttN,rN]
% % ssmPar      - parameter of the state space model
% % Np          - number of particles
% % varObs      - variance of observation noise
% Output
%  %  X     - samples of paths of state variable
%  %  ess   - effectitve sample size of the PF traj in each MCMC step
%  %  theta - samples of theta
% 
% Reference
%   F. Lindsten and M. I. Jordan T. B. Schön, "Ancestor sampling for
%   Particle Gibbs", Proceedings of the 2012 Conference on Neural
%   Information Processing Systems (NIPS), Lake Taho, USA, 2012.
%
% Last updated by Fei Lu, 2019-10-31 (set plotON for demonstration)


X     = zeros(ssmPar.d,tN,numMCMC);
K     = length(theta0); 
theta = zeros(K,numMCMC);

ttN   = length(obs);     % number of SMC time steps = obs Times
ess   = zeros(1,numMCMC);
ssmPar.bounds = prior.statebounds; % lower and upper bounds of state

% Initialize the state by running a PF
[particles,w] = cpf_as(obs,t0t1, ssmPar,dt,Np,varObs,theta0, X(:,:,1),semiEM,1);
ess(1)        = 1/(w(:,end)'*w(:,end));  % effective sample size

% Draw J
J = find(rand(1) <= cumsum(w(:,ttN)),1,'first');
X(:,:,1) = particles(:,:,J);  
if sp ==0 % estimate parameter+ state; Otherwise, estimate state with TruePar  
    theta(:,1) = sampleTheta_Bayes(X(:,:,1),ssmPar,K,dt,prior); 
    theta0     = theta(:,1)';  
end

  if exist('ProgressOn', 'var') == 0; ProgressOn=0; end
      
% Run MCMC loop
reverseStr = []; plotON = 0; 
for k = 2:numMCMC
   if ProgressOn ==1   
    reverseStr = displayprogress(100*k/numMCMC, reverseStr);
   end
    % Run CPF-AS
    if k>3; plotON=0; end
    [particles,w] = cpf_as(obs,t0t1,ssmPar,dt,Np,varObs,theta0,X(:,:,k-1),semiEM,plotON);
    ess(k)    = 1/(w(:,end)'*w(:,end));  % effective sample size
    % Draw J (extract a particle trajectory)
    J = find(rand(1) < cumsum(w(:,ttN)),1,'first');
    X(:,:,k) = particles(:,:,J);  
    
    if sp ==0  % estimate parameter+ state; Otherwise, estimate state with TruePar
        theta0   = sampleTheta_Bayes(X(:,:,k),ssmPar,K,dt, prior);
        theta(:,k)= theta0;
        theta0    = theta0';
    end
end

end


%-------------------------------------------------------------------
function reverseStr = displayprogress(perc,reverseStr)
msg = sprintf('%3.1f', perc);
fprintf([reverseStr, msg, '%%']);
reverseStr = repmat(sprintf('\b'), 1, length(msg)+1);
end


