function [obs, ssmPar,Utrue] = generateData(prior,ssmPar,Datafilename, saveON, plotON,semiEM)
% generate simulated data  
% The state space model: 
%   x_{t+1} = x_t + dt*g(\theta,x_t)+ v_t,
%   y_t     = x_t + W_t,
% where 
%    g(theta,x) = theta* terms 
% is provided by nl_fn(theta,u). 
% 
% Last updated: 2018/7/21

% thetaTrue = [0.1,-1,0.4,-1,0.2]; % original test
if prior.flag ==1
    temp = (prior.mu)';    % should sample from prior later
    thetaTrue = [exp(temp(1)), temp(2:4), -exp(temp(5))]; % 1x5 vector
elseif prior.flag ==0     % Gaussian
    temp = prior.mu + randn(size(prior.mu)).*sqrt(prior.sigma); thetaTrue = temp'; 
elseif prior.flag ==2     % uniform
    temp = prior.mu +randn(size(prior.mu)).*(prior.ub-prior.lb); thetaTrue = temp'; 
end
tN     = ssmPar.tN; 
dt     = ssmPar.dt;
stdObs = ssmPar.stdObs; % std of observation noise % ======== NILS: stochastic forcing =================
% is on the order sqrt(dt)*stdF ~ 0.01; if StdObs is one order higher, it seems impossible
% for me to infer the true process, as the observation noise cancels out the stochastic 
% forcing signal in the process
  
stdF = ssmPar.stdF; 
 
%% forward simulation
d      = ssmPar.d;      % dimension of state
U0     = 0*randn(d,1) + ones(d,1); 
if exist('semiEM','var') ==0; semiEM =0; end
Utrue  = forward_solver(thetaTrue,dt,tN,U0,stdF,semiEM);

%% generate noise observation at each time
 t0t1= [(1:tN)',(2:tN+1)']; 
% restrictions on time intervals: ordered, no overlap,
ttN = length(t0t1(:,1));   

F   = fnObs(Utrue, t0t1);  
obs = F + stdObs* randn(size(F)); 

if nargin>4 && plotON==1
    tt = dt*(2:tN+1); 
    figure; plot(tt,Utrue(2:end),'b'); hold on; 
            plot(tt,obs,'ko','linewidth',1); hold off
            xlabel('Time'); ylabel('State'); title('True and observed signal'); 
            legend('true','Obs'); 
end

% update ssmPar 
ssmPar.t0t1       = t0t1; 
ssmPar.thetatrue  = thetaTrue;  % a row vector

if saveON == 1
save(Datafilename,'obs','ssmPar','prior','thetaTrue','dt','Utrue','t0t1','stdObs','semiEM'); 
end

return

