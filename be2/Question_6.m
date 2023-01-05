%
%          Y. Diouane
%          BE: Problemes Inverses
%          ISAE-SUPAERO
%
%  Question 6 
%
clear all; close all;
global nt dt dx zs source p_obs cmin cmax


linewidth = 1.0;
fontsize = 14;

nt = 1500;
%To avoid numerical dispersion and instability, dx < cmin/(5*fpeak) and dt < 0.6*dx/cmax.
dx = 5;
dt = 0.001;
nz = 301;
c_true = 2000*ones(nz,1);
c_true(51:80) = 1600;
c_true(121:150) = 2500;
cmin=min(c_true);
cmax=max(c_true);
zs = 5; % The source position
list_f = [2,5,10]
t = (0:nt-1)'*dt;
z = (0:nz-1)'*dx;
print=0;
c = 2000*ones(nz,1);
fprintf('\n**********************\nRunning: Multiscale Approach. \n');
n_iter=30;
tolg=1e-10; % Tol. of the gradient.
for j=1:length(list_f) % For each Level
    f1=list_f(j);
    source=SourceTerm(nt,f1,dt,print);
    p0 = zeros(nz,1);
    p1 = p0;
    p0(zs) = source(1);
    p1(zs) = source(2);
    [p_true,p_all]=ForwardProblem(nz,nt,source,c_true,zs,dx,dt,p0,p1,print);
    sig = 0.1;
    p_obs=p_true+sig*randn(nt, 1); % possiblly with noise
  
    fprintf('\n**********************\nRunning Algo Descent : Level with the peak frequency f= %d\n',f1);
    fprintf('it\t obj\t\t norm(d)\t step\n');
    for i=1:n_iter 
        [objfun,g]=CostFunc_FWI(c);
        alpha=100/norm(g);
        d = -g;
        c=  c + alpha*d;
        fprintf('%d\t %.3e\t %.3e\t %.3e \n',i,objfun,norm(g), alpha);
        figure(3)
        if(j==1)
            plot(z,c_true,'b',z,c,'k'); % f=2
        elseif(j==2)
            plot(z,c_true,'b',z,c,'c'); % f=5
        else
            plot(z,c_true,'b',z,c,'r'); % f=10
        end
        legend('True Velocity','Inverted Velocity');
        xlabel('Depth (m)');ylabel('Velocity (m/s)');
        title(['Peak Frequency f= ', num2str(f1)]);
        axis([0 1500 1000 3000]);
        drawnow;
   
        if(norm(g) < tolg) % Define your stopping criteria
            break;
        end
    end
    if(j==1)
        c_f2=c;
    end
    if(j==2)
        c_f5=c;
    end
    if(j==3)
        c_f10=c;
    end
end
%
% Diplay all inverted velocities at once
%
figure(4)
clear print
plot(z,c_true,'b',z,c_f2,'k',z,c_f5,'c',z,c_f10,'r', 'LineWidth', 1.4);
grid minor
leg = legend('True Velocity','Inverted Velocity (f=2)','Inverted Velocity (f=5)','Inverted Velocity (f=10)');
xl=xlabel('Depth (m)');
yl=ylabel('Velocity (m/s)');
set(yl, 'Interpreter', 'latex', 'fontsize', fontsize , 'LineWidth', linewidth);
set(xl, 'Interpreter', 'latex', 'fontsize', fontsize , 'LineWidth', linewidth);
set(leg, 'Interpreter', 'latex', 'fontsize', fontsize , 'LineWidth', linewidth, 'Location', 'southeast', 'NumColumns', 1); % southeast bestoutside

%title('Inversion using  the Multiscale Approach');
axis([0 1500 1400 2600]);
clear print
titles = ['q6_nf:', num2str(length(list_f))];
print(4, titles, '-dpng')
