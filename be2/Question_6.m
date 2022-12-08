%
%          Y. Diouane
%          BE: Problemes Inverses
%          ISAE-SUPAERO
%
%  Question 6 
%
clear all; close all;
global nt dt dx zs source p_obs cmin cmax

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
list_f=[2,5,10]
t = (0:nt-1)'*dt;
z = (0:nz-1)'*dx;
print=0;
c = 2000*ones(nz,1);
fprintf('\n**********************\nRunning: Multiscale Approach. \n');
n_iter=30;
tolg=1e-10; % Tol. of the gradient.
for j=1:length(list_f) % For each Level
    f1=list_f(j); % Choose a peak frequency
    %
    % Compute p_obs associated with f1
    % (To do)
    p_obs= 
    fprintf('\n**********************\nRunning Algo Descent : Level with the peak frequency f= %d\n',f1);
    fprintf('it\t obj\t\t norm(d)\t step\n');
    for i=1:n_iter % Launch optimization
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
        % PRINT
        %
        fprintf('%d\t %.3e\t %.3e\t %.3e \n',i,objfun,norm(g), alpha);
        %
        % DISPLAY FIGURE
        %
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
        %
        % Check stopping criteria
        %
        if( ) % Define your stopping criteria
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
plot(z,c_true,'b',z,c_f2,'k',z,c_f5,'c',z,c_f10,'r');
legend('True Velocity','Inverted Velocity (f=2)','Inverted Velocity (f=5)','Inverted Velocity (f=10)');
xlabel('Depth (m)');ylabel('Velocity (m/s)');
title('Inversion using  the Multiscale Approach');
axis([0 1500 1000 3000]);