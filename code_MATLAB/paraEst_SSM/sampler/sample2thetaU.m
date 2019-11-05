function [thetaS,Usample]=sample2thetaU(sample,Nb,tN,thetaL)
   % get back to theta and U from sample
       [~,ns] = size(sample); 
       thetaS = sample(1:thetaL,:);  % samples of theta
       Usample = zeros(Nb,tN,ns);
       
       for i=1:ns
            temp = sample(thetaL+1:end,i);
            Usample(:,:,i) = reshape(temp,[Nb,tN]);
       end    
   end