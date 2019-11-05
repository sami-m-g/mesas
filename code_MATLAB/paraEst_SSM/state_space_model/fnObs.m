function  F = fnObs(U, t0t1)
% Observation function 
%    1. Time intervals ordered, non-overlap
% The function is an integral in region and from time t0 to t1 
%       fn = \sum_{t=t0+1}^{t1} U(t) /(t1-t0)
% Input
%    U       - the coefficients in the FEM
%    t0t1    - start and end time of each region       tNx2
% Output
%    F     - the value of the function
%  
% Last updated by Fei Lu, 2018/7/21

tN  = length(t0t1(:,1)); 
F   = zeros(tN,1);
for tt =1:tN
    tindx    = t0t1(tt,1)+1:1:t0t1(tt,2);  
    utemp    = U(tindx); 
    F(tt)    = sum(utemp)/length(utemp);
end

return