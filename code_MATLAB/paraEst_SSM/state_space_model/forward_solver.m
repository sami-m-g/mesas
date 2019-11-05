function [U,U1det] = forward_solver(theta,dt,tN,U0,stdF,semiEM)
% 1D equation forward by Euler scheme
% The state space model: 
%   x_{t+1} = x_t + dt*g(\theta,x_t)+ W_t,
%   y_t     = x_t + V_t,
% By Fei Lu, Nils Weitzel, 18/7/21

Nnodes = 1;  
U      = zeros(Nnodes,tN+1);   
%% Time integration  
U(:,1) = U0;     % Initial Condition
if stdF>0         % random increment 
    Winc = stdF* randn(Nnodes,tN+1);
else 
    Winc = zeros(Nnodes,tN+1); 
end

periodfn =@(x) 0.1*sin(x*2*pi); 

if semiEM ==0
    % the nonlinear term
    % theta = [0,0,1,1,-0.1];  % parameter p0,...,p4
    nlfn = @(u) nl_fn(u,theta);
    for n = 2:tN+1
        % Non-linear forces: at the gravity center of the triangle
        nl      = nlfn(U0);  % nolinear term
        nl      = nl+ periodfn(n*dt); 
        % Stochastic forcing:
        stocf = sqrt(dt)*Winc(n);
        
        if n==2; U1det = U0 + dt* nl; end
        U1     = U0 + dt* nl + stocf;
        U(:,n) = U1;    U0     = U1;
        if abs(U1)>1e8; fprintf('\n Solution blows up at time %d \n',n); return; end
    end
elseif semiEM ==1
    nlfn = @(u) nl_fn(u,theta); 
    cdt  = 1- theta(2)*dt; 
    for n = 2:tN+1
        nl      = nlfn(U0);         % nolinear term        
        nl      = nl+ periodfn(n*dt); 
        stocf = sqrt(dt)*Winc(n);   % Stochastic forcing:      
        if n==2; U1det = U0 + dt* nl; end
        U1     = U0 + (dt* nl + stocf)/cdt;
        U(:,n) = U1;    U0     = U1;
        if abs(U1)>1e3; fprintf('\n Solution blows up at time %d \n',n); return; end
    end
end

return



