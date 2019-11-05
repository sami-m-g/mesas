function [nl,terms] = nl_fn(u,theta,t)
% nonlinear function in the stochastic energy balance model
% Input
%   theta    - coefficients of the nonlinear terms
%   u        - solution u 
%   t        - time
% Output
%   nl       - value of the nonlinear term at time t

numtheta = length(theta); 

switch  numtheta
    case 5    %%%---- 5 terms
        terms = [ones(size(u)); u; u.^2; u.^3; u.^4];
        nl = theta*terms; 
    case 4    %%% ---- 4 terms
        terms = [ -(u-1); u.*(u-1);u.^2.*(u-1); u.^3.*(u-1)];  %  [(-1; u; u.^2; u.^3),*(u-1); -tan]
        % terms = [ -(u-1); exp(u).*(u-1); 10*cos(2*pi*u).*(u-1); u.^3.*(u-1)];
        nl = theta*terms - 10^(-4)*tan(2*pi/2*(u-1)); 
    case 3    %%%  3 terms
        
        %%%--- % Fei's function
%        terms = [ -(u-1);  -1*cos(4*pi*u).*(u-1); u.^3.*(u-1)];
%        nl = theta*terms - 10^(-4)*tan(2*pi/2*(u-1)); 
        
        %%%---  Adam's function
        terms = [ones(size(u)); u; u.^4];   
        nl    = theta*terms; 
        % nl    = nl - 10^(-4)*tan(2*pi/2*(u-1));
    case 2
        terms = [ones(size(u)); u.^1];   
        nl    = theta*terms; 
end



return