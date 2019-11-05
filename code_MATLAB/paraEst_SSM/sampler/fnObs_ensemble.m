function  [F] = fnObs_ensemble(UU)
% evaluate the observation function at an ensemble at ONE time interval t0:t1
%   -- used in SMC for the computation of weight/target distribution
% The observation function is an integral in region and from time t0 to t1 
%       fn = mean (UU,2)
% Input
%    U       - an ensemble of the solution in [t0,t1]:  size [Dx,Nt,Np]
% Output
%    F     - the value of the function,   1 x Np
%       
F  = mean(UU,2); 
F  = reshape(F,[],size(UU,3)); 

return