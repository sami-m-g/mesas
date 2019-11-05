function plot_g(thetatrue,lowupbounds)
% plot the nonlinear function g(theta,x)
 a = lowupbounds(1); b = lowupbounds(2); dx = (b-a)/200; 
x = a:dx:b; 
gvalue = nl_fn(x,thetatrue); 
plot(x,gvalue,'linewidth',1); hold on; 
% plot(x,0*x,'-.');  
xlabel('x'); ylabel('g(x)'); title('The nonlinear function g(\theta,x)' );
% axis([0,2,-600,600]);
print('fig_new_g','-depsc');

end
