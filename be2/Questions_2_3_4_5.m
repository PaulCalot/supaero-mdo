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
values = [1950 2000 2050 2100];
[p_true,p_all]=ForwardProblem(nz,nt,source,c_true,zs,dx,dt,p0,p1,print);



lbl = cell(1,length(values)+1);
t = (0:nt-1)'*dt;
z = (0:nz-1)'*dx;
figure(4)
hold on
plot(z,c_true,'b', 'LineWidth', 1.5)
lbl{1} = 'True Velocity';
quant_lab = 'c0';
%
% For this starting point c_0
%
for v = 1:length(values)
    value = values(v);
    
    sig = 0.1;
    p_obs=p_true+sig*randn(1)*p_true; % possiblly with noise

    c = value*ones(nz,1); 

    fprintf('\n**********************\nRunning Algo Descent : peak frequency f= %d\n',f1);
    fprintf('it\t obj\t\t norm(d)\t step\n');
    n_iter = 200;
    tolg=1e-10; % Tol. of the gradient.
for i=1:n_iter 
    %
    %    evaluate the objective function and its gradient
    %    (To do)
    [objfun,g]= CostFunc_FWI(c);
        %
        %  Update the step size alpha and the velocity profile c
        %  (To do ) 
        alpha = 100/norm(g);
        %
        % Compute the search direction 
        %
        d = -g;
        %
        % Update the velocity profile
        %
        c = c + alpha*d;
        %
        %
        fprintf('%d\t %.3e\t %.3e\t %.3e \n',i,objfun,norm(g), alpha);
        %
        % DISPLAY FIGURE
        %
%     figure(3)
%     plot(z,c_true,'b',z,c,'r');legend('True Velocity','Inverted Velocity');
%     xlabel('Depth (m)');ylabel('Velocity (m/s)');
%     title('Low frequency inversion');
%     axis([0 1500 1000 3000]);
%     drawnow;
    %
    % Check stopping criteria 
    %
    if(norm(g) < tolg) % Define your stopping criteria
        break;
    end
end
plot(z,c,'LineWidth', 1.1);
lbl{v+1} = ['$c_0 $ = ', num2str(value)];

axis([0 1500 1000 3000]);
end
% ttl = title('Inversed Velocity using  Descent Algo');
% set(ttl, 'Interpreter', 'latex', 'fontsize', 20);
grid minor
xl = xlabel('Depth (m)');
yl = ylabel('Velocity (m/s)');
set(yl, 'Interpreter', 'latex', 'fontsize', 16 , 'LineWidth', 1.5);
set(xl, 'Interpreter', 'latex', 'fontsize', 16 , 'LineWidth', 1.5);
leg = legend(lbl);
set(leg, 'Interpreter', 'latex', 'fontsize', 14 , 'LineWidth', 1.5, 'Location', 'southeast');
clear print
titles = ['Quest4 for var', quant_lab];
print(4, titles, '-deps')

%figure(5)
%plot(residual);title('Residual over iterations');
%xlabel('Iteration');ylabel('Residual');
