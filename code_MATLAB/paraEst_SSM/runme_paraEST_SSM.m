% sample the posterior of the states and the parameters
%  --- Method: particle Gibbs with Ancestor Sampling (PGAS)
%  --- only considers parameters in the nonlniear term;
%      the variance of the noise in the SEB is assumed known
%
% Last updated: Fei Lu, 2019-11-1
% NEXT steps:
%  - Gaussian process regression code 
%  - formulation: see if the problem can be formulated into linear in para;
%     >>> needed for either Particle Filter or Gaussian Process regression
  

clear all; close all; 
restoredefaultpath;  addpaths; 

%% settings: state model information
stdObs = 0.10;    noise = '_n01'; 
dt   = 0.01;    % unit year
tN   = 100;    % number of time steps
stdF = 0.1;     

ssmPar.stdF       = stdF;         % std of stochastic force
ssmPar.stdObs     = stdObs;       % std of observation noise
ssmPar.d          = 1;            % dimension of state
ssmPar.dt         = dt;           % time step size
ssmPar.tN         = tN;           % number of time steps 

settings = [', tN = ', num2str(tN), ', stdObs = ', num2str(stdObs)];  

%% Settings: prior of parameters 
prior.flag  = 0; nldeg = '014';
if prior.flag == 2;      type = 'Unif';   % 0 Gaussian, 2 uniform
elseif prior.flag ==0;   type = 'Gauss';   % 0 Gaussian, 2 uniform
end   
sp = 0;   % 0 estimate the parameter; 1, estimate states with true para
if sp ==1; type1 = strcat('s_',type);  % estimate state, uring TRUE paramter      
else;      type1 = strcat('sp_',type); % estimate state and paramter  
end

path = ['/outputs/newgDeg',nldeg,'/data/'];  % path to save data
fprintf(['Running: Prior= ', type, ', terms= ', nldeg, settings, '\n \n']); 

figname = ['_tN',num2str(tN),type1,noise]; 
sampleFilename = strcat(pwd,path,'SampleNewIC',figname,'.mat'); 
Obsdatafile = strcat(pwd,path,'ObsData',figname,'.mat');   
 
%%  Generate observation data and Sample the posterior by particle MCMC: PGAS
%% generate Observation data
if exist(Obsdatafile,'file') == 0
    % set prior 
    saveON = 1;
    theta_alpha;     % load thetatrue, lower and upperbound, thetaStd  ***********
    thetaTrue = thetatrue;
    prior.statebounds=[0.6,1.4];  % lower and upper bounds of state
    if  prior.flag == 0       % Gaussian prior
        prior.mu  = thetaTrue';   prior.sigma =1*thetaStd'.^2;  % sigma_i^2
    elseif prior.flag ==2     %%% === Uniform prior
        prior.mu  = thetaTrue';  % Here prior.mu used to pass the true value.
        prior.lb  = thetaLo';     prior.ub  = thetaUp';
    end
    
    plotON =1;  semiEM = 1;
     generateData(prior,ssmPar,Obsdatafile,saveON,plotON,semiEM);
end
load(Obsdatafile);

%%% plot the histogram of Utrue and estimated nonlinear fcn g
figure; histogram(Utrue,100); print([pwd,path,'hist',figname],'-depsc','-r600'); title('Histogram: state variable');

%% MLE from true and obs data
[mleTrue,cov1,sigmaMLE]      = MLE_truestate(Utrue,ssmPar,semiEM);
[mleObs,mleMat,sigmaObs] = MLE_truestate(obs',ssmPar,semiEM);
if     length(mleObs)==2; type = ' %2.4f %2.4f ';
elseif length(mleObs)==3; type = ' %2.4f %2.4f %2.4f ';
elseif length(mleObs)==4; type = ' %2.4f %2.4f %2.4f %2.4f  ';
end
svdmatrix = '   svd matrix:'; sigmaM = 'sigma  %2.4f ';
fprintf([' True      :', type, '\n'], ssmPar.thetatrue);
fprintf(['MLE(Utrue):', type, svdmatrix, type, sigmaM,'\n'], mleTrue', svd(cov1),sigmaMLE);
fprintf(['MLE(obs)  :', type, svdmatrix, type, sigmaM,'\n'], mleObs',svd(mleMat),sigmaObs);

keyboard;     
%%  %%%%%  Sample the posterior by particle MCMC: PGAS
 if exist(sampleFilename,'file') == 0  % if no sample file, generate samples   
    Np      = 5;         % number of particles in SMC
    numMCMC = 1000;         % length of MCMC chain
    burnin  = 0.3*numMCMC;  % burnin 
    theta0  = sp*thetaTrue;   % inital guess
    varObs  = stdObs^2; 
    % Run the algorithms:   size(Usample)=[ Dx,tN,numMCMC]; theta: [K,numMCMC]
    fprintf('Running PGAS (N=%i). Progress: ',Np); tic;
    [Usample,ess,theta] = pgas_statePar(numMCMC,obs,t0t1,ssmPar,Np,dt,...
        tN+1,theta0,varObs, prior,sp,semiEM,1);   % size(theta)=K x numMCMC 
    timeelapsed = toc;
    fprintf(' Elapsed time: %2.2f sec.\n',timeelapsed); 
 
    clearvars -except timeelapsed  ess theta  numMCMC Np dt burnin ssmPar ...
                      prior mleObs mleTrue mleMat Obsdatafile Usample type ...
                      thetaTrue   Utrue   sampleFilename figname path ;               
    save(sampleFilename);
end

%% represent the results
load(sampleFilename);load(Obsdatafile);
    fprintf(['   True parameters:', type,'\n'], ssmPar.thetatrue); 
    fprintf(['   MLE from Obs:   ', type,'\n'], mleObs');
    fprintf(['   Posterior mean: ', type,'\n'], mean(theta,2)');
    fprintf(['   Prior     mean: ', type,'\n'], prior.mu');
    covPost = cov(theta'); fprintf('Cov Posterior\n'); disp(covPost); 
    fprintf('Cov Prior\n'); disp(diag(prior.sigma));
    

%%% plot estimated nonlinear fcn g
mUobs = mean(obs); stdUobs = std(obs); bds = [mUobs-3*stdUobs, mUobs+3*stdUobs];
figure; plot_g(thetaTrue,bds); hold on; 
plot_g(mleTrue',bds); plot_g(mleObs',bds); plot_g(mean(theta,2)',bds); 
legend('True','MLE','MLEobs','Posterior'); print([pwd,path,'fcn_g',figname],'-depsc'); title('Function g and estimator');


% coverage frequency: good for the highD state
quantile = 0.90;   
covfreq  = coverageFrequency(sampleFilename,quantile,figname,path);


%% plot the marginal of the posterior, and trajectories: of para and state
[sMean,sStd,theta]= Plots_marginal(sampleFilename,figname,path,Obsdatafile);






