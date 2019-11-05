function plotdensity(sampleArray, trueV,titl,prior)
% plot the marginal densities given by samples

K = length(trueV);
for kk = 1:K
  subplot(K,1,kk); 
  samples = sampleArray(kk,:); 
  smean   = mean(samples); % std1 = std(samples); 
  [f,xi]  = ksdensity(samples);
  plot(xi,f); hold on; 
  tmp=get(gca,'ylim');
  plot([trueV(kk) trueV(kk)],tmp*1.5,'k-*');
  plot([smean smean],tmp*1.5,'k-d'); 
  if exist('prior','var')
      plot([prior.mu(kk) prior.mu(kk)],tmp*0.8,'k-x');
      lb = -sqrt(prior.sigma(kk)) +prior.mu(kk);
      ub = sqrt(prior.sigma(kk))+ prior.mu(kk);
      plot([lb,lb],tmp*0.5,'k-x'); plot([ub,ub],tmp*0.5,'k-x');
  end
end

 subplot(K,1,1);  title(titl,'FontSize',16);
 legend({'posterior','true value','posterior mean','prior mean/std'},'FontSize',16)