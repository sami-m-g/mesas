
%--------------------------------------------------------------------------
function [x,w] = cpf_as(obs,t0t1,ssmPar,dt,Np,varObs,theta,X1,semiEM,plotON)
% Conditional particle filter with ancestor sampling
% Input:
%   obs  - measurements: size = [# of time intervals, # regions]
%                            = [length(t0t1(:,1)), length(regions(1,:))]
%   varObs- measurement noise variance
%   Np    - number of particles
%   X1     - conditioned particles - if not provided, un unconditional PF is run

stdF     = ssmPar.stdF; stdF2 = stdF.^2; 

if exist('semiEM','var') ==0; semiEM =0; end
forward  = @(U0,dt,tn) forward_solver(theta,dt,tn,U0,stdF,semiEM); 

conditioning = 1; % (nargin > 9);
[Dx,T,~]     = size(X1); 
X_Np         = reshape(X1,[Dx,T]);  % size(X1)=[ Dx,T,1]->>[Dx T]. Ref path
x        = zeros(Dx,T, Np); % Particles

x(:,1,:) = 1+0.05*randn(Dx,Np);    % positive initial condition

x(:,:,Np) = X1(:,:,1); % Set the Np:th particle according to the conditioning

ttN = length(obs);   % number of SMC time steps = obs Times 

a  = zeros(Np,ttN);   % Ancestor indices
w  = zeros(Np, ttN);  % Weights

tt = 1;      % observation time index 

% NOTE on t0t1: here we assume t0t1, an Lx2 array, are start and end times 
% of time intervals, and Importantly, each start = the previous end, i.e. 
%            t0t1(l,1) = t0t1(l-1,2) 
while tt <= ttN
    t0= t0t1(tt,1); t1 = t0t1(tt,2);  steps =t1 - t0; 
    if tt==1   %%% initialization up to the first t1 
        xt0 = reshape(x(:,t0,:), Dx,Np); ind = 1:Np; 
        [xpred,~] = forward_ensemble(forward,xt0,dt,steps,ind); % size(xpred)=[Dx,steps+1, Np]
        x(:,t0+1:t1,:) = xpred(:,2:end,:); 
    else % % tt>1 
        ind = myResampling(w(:,tt-1));
        ind = ind(randperm(Np));        
        xt0 = reshape(x(:,t0,:), Dx,Np);  % resampled particles as inital condition
        [xpred,u1det] = forward_ensemble(forward,xt0,dt,steps,ind); % -Checked: 6/19-
                     % size(u1det) = [Dx,Np]
                     % u1det is the 1-step deterministic forward of the particles
        x(:,t0+1:t1,:) = xpred(:,2:end,:);   % NEW: removed t0 overlap 
        if(conditioning)
            x(:,t0+1:t1,Np) = X_Np(:,t0+1:t1); % Set the Np-th particle as conditioning
            % Ancestor sampling
            X_Np1= X_Np(:,t0+1);    % NEW: only one-step probability needed
            wexp = multiPriorWt(X_Np1,u1det,dt,stdF2); 
            wexp = wexp + logweights'; % use the exponents to avoid singular
            w_as = exp(wexp- max(wexp));  
            
            w_as = w_as/sum(w_as);
            temp = find(rand(1) <= cumsum(w_as),1,'first');
             if isempty(temp); keyboard; end
            ind(Np) = temp;
        end
        % Store the ancestor indices --- starting from tt=2
        a(:,tt) = ind;
    end
    if exist('plotON','var') && plotON ==1
        figure(17); clf; titl='Before traj-resample'; plotTrajpf(x,1,titl); pause(0.05); 
    end
    % Compute importance weights:
    UU     = x(:,t0:t1,:); % Let us Right-point approximation for now
    ypred  = fnObs_ensemble(UU); % size=[rN,Np]
    logweights = -1/(2*varObs)*sum( (obs(tt,:)' - ypred).^2,1); 
    const = max(logweights); % Subtract the maximum value for numerical stability
    weights = exp(logweights-const);
    w(:,tt) = weights/sum(weights); % Save the normalized weights
    tt = tt+1; 
end

% Generate the trajectories from ancestor indices
ind = a(:,ttN);
for tt = ttN-1:-1:1
    t0= t0t1(tt,1); t1 = t0t1(tt,2);
    if tt>1
        x(:,t0+1:t1,:) = x(:,t0+1:t1,ind);
    else
        x(:,t0:t1,:) = x(:,t0:t1,ind);
    end
    ind = a(ind,tt);
end
end


%-------------------------------------------------------------------
function  wexp = multiPriorWt(X_Np1,xpred_det,dt,stdF2)  % UPDATED
% compute likelihood weight in ancester sampling  >> only one step needed
%   X_Np1    - reference trajectory at t0+1 (Corrected from t0:t1) size = [Dx,1]
%   xpred_det- predicted by deterministic forward: size = [Dx,Np]
%   chol_PAB  - Cholesky facotorization of the covariance matrix
% Formular
%      U_{n+1} = pred_det + stdF* N(0, dt M) 

[~,Np] = size(xpred_det);
wexp     = zeros(Np,1);
for nn = 1:Np
    diff  = X_Np1 - xpred_det(:,nn);  
    pdiff = 1 * diff;   
    wexp(nn) = - sum(pdiff.^2)/(2*stdF2*dt); 
end
end
