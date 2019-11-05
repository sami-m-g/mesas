function [X,ess] = pgas_state(numMCMC,obs,t0t1,regions,femPar,Np,dt,tN,theta,r0)
% Runs the PGAS algorithm,
%  % numMCMC     - number steps in MCMC
%  % obs         - observations:   size =[ttN,rN]
%  % t0t1,regions- cooresponding to obs
%  % x0          - initial condition
%  % Np          - number of particles
%  % q0          - covariance matrix of noise in forward model
%  % r0          - std of observation noise
% Output
%  % samples of X
% Reference
%   F. Lindsten and M. I. Jordan T. B. Schön, "Ancestor sampling for
%   Particle Gibbs", Proceedings of the 2012 Conference on Neural
%   Information Processing Systems (NIPS), Lake Taho, USA, 2012.
%
% The function returns the sample paths of (x_{1:T}).

Prec_chol = femPar.Prec_chol; 
Dx  = length(Prec_chol(1,:)); 
X   = zeros(Dx,tN,numMCMC);


ttN = length(t0t1(:,1));   % number of SMC time steps = obs Times
ess = zeros(1,numMCMC);

% Initialize the state by running a PF
[particles,w,sFlag] = cpf_as(obs,t0t1,regions,femPar,dt,Np,r0, X(:,:,1));
ess(1)    = 1/(w(:,end)'*w(:,end));  % effective sample size
% Draw J
J = find(rand(1) <= cumsum(w(:,ttN)),1,'first');
X(:,:,1) = particles(:,:,J); 


% Run MCMC loop
reverseStr = [];
for(k = 2:numMCMC)
    reverseStr = displayprogress(100*k/numMCMC, reverseStr);
    
    % Run CPF-AS
    [particles,w,sFlag] = cpf_as(obs,t0t1,regions,femPar,dt,Np,r0, X(:,:,k-1));
    ess(k)    = 1/(w(:,end)'*w(:,end));  % effective sample size
    % Draw J (extract a particle trajectory)
    J = find(rand(1) < cumsum(w(:,ttN)),1,'first');   
    X(:,:,k) = particles(:,:,J);    
end

end




%--------------------------------------------------------------------------
function [x,w,sFlag] = cpf_as(obs,t0t1,regions,femPar,dt,Np,r0,X1)
% Conditional particle filter with ancestor sampling
% Input:
%   obs  - measurements: size = [# of time intervals, # regions]
%                            = [length(t0t1(:,1)), length(regions(1,:))]
%   t0t1 - start and end of time intervals  ttNx2; connected (for SMC)
% regions- nodes in regions, size= [ttN*rN, eN]
%   r0    - measurement noise variance
%   Np    - number of particles
%   X     - conditioned particles - if not provided, un unconditional PF is run

Prec_chol = femPar.Prec_chol; 
elements3 = femPar.elements3; 
TkVol     = femPar.tri_areas/6; 
forward   = femPar.forward; 

conditioning = (nargin > 8);
[Dx,T,~]     = size(X1); 
X_Np         = reshape(X1,[Dx,T]);   %  size(X1) is [ Dx T, 1]--->> [Dx T]
x  = randn(Dx,T,Np); % Particles
x(:,:,1) = randn(Dx,T);    % Deterministic initial condition
x(:,:,Np) = X1(:,:,1); % Set the Np:th particle according to the conditioning


ttN = length(t0t1(:,1));   % number of SMC time steps = obs Times 
rN  = length(regions(:,1)) / ttN; 

a  = zeros(Np, ttN);   % Ancestor indices
w = zeros(Np, ttN);    % Weights


tt=1;      % observation time index 
sFlag = 0 ;  % singular indicator
singN    = 0; % singular increment distribution
while tt <= ttN
    t0= t0t1(tt,1); t1 = t0t1(tt,2);  steps =t1 - t0; 
    if (tt~=1)   
        ind = myResampling(w(:,tt-1));
        ind = ind(randperm(Np));
 %        x   = x(:,:,ind);               % resampling the whole trajectory
        xt0 = reshape(x(:,t0,:), Dx,Np);
        xpred = forward_ensemble(forward, xt0,dt,steps);    % size(xpred) = [Dx,Np]
        x(:,t0:t1,:) = xpred(:,:,ind);
        if(conditioning)
            x(:,t0:t1,Np) = X_Np(:,t0:t1); % Set the Np:th particle according to the conditioning
            % Ancestor sampling
            X_Np1= X_Np(:,t0:t1);
            wexp = multiPriorWt(X_Np1,xpred,Prec_chol,dt); %multistep prior weight
            wexp = wexp - max(wexp); 
            m = exp( wexp/2);
            w_as = w(:,t-1).*m;      % m and w(:,t-1) can be singular!!!!
            if sum(w_as) ==0         % count tries of resample
                singN = singN +1;    % disp(singN);
                if singN >10^3
                    disp('Singular increment distribution. Stopped.'); 
                    sFlag =1;       return;     end
                continue; 
            end 
            w_as = w_as/sum(w_as);
            temp = find(rand(1) <= cumsum(w_as),1,'first');
            if isempty(temp); keyboard; end
            ind(N) = temp;
        end
        % Store the ancestor indices
        a(:,tt) = ind;
    end
    % Compute importance weights:
    UU     = x(:,t0:t1,:);
    region1= regions(rN*(tt-1)+(1:rN),:);   
    ypred = fnObs_ensemble(UU, region1, elements3, TkVol); % size=[rN,Np]
    logweights = -1/(2*r0)*sum( (obs(tt,:)' - ypred).^2,1); 
    const = max(logweights); % Subtract the maximum value for numerical stability
    weights = exp(logweights-const);
    w(:,tt) = weights/sum(weights); % Save the normalized weights
    tt = tt+1; 
end

% Generate the trajectories from ancestor indices
ind = a(:,ttN);
for(tt = ttN-1:-1:1)
    x(:,tt,:) = x(:,tt,ind);
    ind = a(ind,tt);
end
end

%-------------------------------------------------------------------
function  wexp = multiPriorWt(X_Np1,xpred,Prec_chol,dt)
% compute multistep prior prior weight in ancester sampling
%   X_Np1    - the reference trajectory t0:t1
%   xpred    - predicted ensemble trajectories
%   Prec_chol- chol of precision matrix; 

[Dx,ttN,Np] = size(xpred);
wexp     = zeros(Np,1);
for nn = 1:Np
    diff  = X_Np1(:,:) - reshape(xpred(:,:,nn),[Dx,ttN]);  
    pdiff = Prec_chol * diff; 
    wexp(nn) = sum(pdiff.^2,1)/dt; 
end

end

%-------------------------------------------------------------------
function reverseStr = displayprogress(perc,reverseStr)
msg = sprintf('%3.1f', perc);
fprintf([reverseStr, msg, '%%']);
reverseStr = repmat(sprintf('\b'), 1, length(msg)+1);
end
