
function  [Bn,Gn] = BnGn(dt,U0,U1,K)
% the terms from state model for likelihood
%           x_{t+1} = x_t + dt*g(\theta,x_t)+ v_t,
%  with g(theta,x) = theta*G(x) 
%  Denoteing the terms  
%        Bn = x(t+1) - x(n)     1x1
%        Gn = dt * G(xn)        Kx1
%  the log likelihood of theta can be written as: 
%        - | Bn -theta*Gn |^2/ sigmaV^2
% Input
%    U0,U1        - U(:,n-1) and U(:,n), or their vector;  NxtN
%    dt           - time step size
%    nlfn         - the nonlinear function  g(u,t)
%    K            - dimension of parameter
% Output
%    Bn Gn        
 
Bn = U1- U0;
[~,terms] = nl_fn(U0,zeros(1,K)); 
Gn = dt*terms; 
end 