%
%
%          Y. Diouane
%          BE: Problemes Inverses
%          ISAE-SUPAERO
%
% For questions 2,3,4, and 5
%
clear all; close all;

global nt dt dx zs source p_obs cmin cmax

nt = 1500;
%
%To avoid numerical dispersion and instability, dx < cmin/(5*fpeak) and dt < 0.6*dx/cmax.
%
dx = 5;
dt = 0.001;
nz = 301;
c_true = 2000*ones(nz,1);
c_true(51:80) = 1600;
c_true(121:150) = 2500;
cmin=min(c_true);
cmax=max(c_true);
zs = 5; % The source position
f1 = 5; % The peak frequency
print=0;
source=SourceTerm(nt,f1,dt,print);
p0 = zeros(nz,1);
p1 = p0;
p0(zs) = source(1);
p1(zs) = source(2);
[p_true,p_all]=ForwardProblem(nz,nt,source,c_true,zs,dx,dt,p0,p1,print);
sig = 0;
p_obs=p_true+sig*randn(1)*p_true; % possiblly with noise
t = (0:nt-1)'*dt;
z = (0:nz-1)'*dx;
%
% For this starting point c_0
%
c = 2000*ones(nz,1); 

fprintf('\n**********************\nRunning Algo Descent : peak frequency f= %d\n',f1);
fprintf('it\t obj\t\t norm(d)\t step\n');
n_iter = 50;
tolg=1e-10; % Tol. of the gradient.
for i=1:n_iter 
    %
    %    evaluate the objective function and its gradient
    %    (To do)
    [objfun,g]=
        %
        %  Update the step size alpha and the velocity profile c
        %  (To do ) 
        alpha=
        %
        % Compute the search direction 
        %
        d = 
        %
        % Update the velocity profile
        %
        c= 
        %
        %
        fprintf('%d\t %.3e\t %.3e\t %.3e \n',i,objfun,norm(g), alpha);
        %
        % DISPLAY FIGURE
        %
    figure(3)
    plot(z,c_true,'b',z,c,'r');legend('True Velocity','Inverted Velocity');
    xlabel('Depth (m)');ylabel('Velocity (m/s)');
    title('Low frequency inversion');
    axis([0 1500 1000 3000]);
    drawnow;
    %
    % Check stopping criteria 
    %
    if( ) % Define your stopping criteria
        break;
    end
end
figure(4)
plot(z,c_true,'b',z,c,'k');legend('True Velocity','Inverted Velocity');
xlabel('Depth (m)');ylabel('Velocity (m/s)');
title('Inversion using  Descent Algo');
axis([0 1500 1000 3000]);
%figure(5)
%plot(residual);title('Residual over iterations');
%xlabel('Iteration');ylabel('Residual');
