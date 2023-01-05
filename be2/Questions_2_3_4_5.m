%
%
%          Y. Diouane
%          BE: Problemes Inverses
%          ISAE-SUPAERO
%
% For questions 2,3,4, and 5
%
clear
close all

global nt dt dx zs source p_obs cmin cmax

analysis_to_make = 2;
noise_on_observed_pressure = 0.2; % eta parameter - 0 for Q2 and 3, 0.1 for Q4 and 5 seems enough
plot_pressure = 0; % 1 to plot pressure
follow_optimisation = 0;
% Analysis 
% 0 : c0
% 1 : sk
% 2 : profile (only first crenel is shifted)
% 3 : smoothing

nt = 1500;
% To avoid numerical dispersion and instability, dx < cmin/(5*fpeak) and dt < 0.6*dx/cmax.
dx = 5;
dt = 0.001;
nz = 301;
zs = 5; % The source position
f1 = 5; % The peak frequency
print=0; % for the forward problem - to print the wave

% Stopping criteria 
n_iter = 500; % max iteration
tolg= 1e-10; % Tol. of the gradient.

t = (0:nt-1)'*dt;
z = (0:nz-1)'*dx;

% figure
figure(4)
hold on
colors = {[0.8500 0.3250 0.0980],... 
[0.9290 0.6940 0.1250],...
[0.4940 0.1840 0.5560],...
[0.4660 0.6740 0.1880],...
[0.3010 0.7450 0.9330],...
[0.6350 0.0780 0.1840]};
% colors = ['k', 'r', 'g', 'm', 'b'];
% colors = ['#0072BD', '#D95319', "#EDB120", "#7E2F8E", "#77AC30"];
linewidth = 1.5;
fontsize = 16;
legend_fontsize = 13;

if(analysis_to_make == 0)
    fprintf('Analysis on c0')
    quant_lab ='c0';
    values = [1950, 1990, 2000, 2010, 2050];
    lbl = cell(1, length(values)+1);
    question_name = 'q2';
elseif(analysis_to_make == 1)
    fprintf('Analysis on sk')
    quant_lab ='sk';
    values = [10, 50, 100, 150, 200];
    lbl = cell(1, length(values)+1);
    question_name = 'q2';
elseif(analysis_to_make == 2)
    fprintf('Analysis on velocity profile')
    quant_lab ='shift';
    values = [-2, 0, 2]; % number of space step - is multiplied by dx - be careful !
    lbl = cell(1, length(values));
    question_name = 'q3';
elseif(analysis_to_make == 3)
    fprintf("Analysis on a priori velocity")
    quant_lab ='smooth';
    values = [0.01, 0.05, 0.1];
    lbl = cell(1, length(values));
    question_name = 'q5';
else
    fprintf('Simple Run')
    quant_lab ='';
    values = 0;
end

if(noise_on_observed_pressure > 0)
    question_name = 'q4';
end

titles = [question_name, '_', quant_lab, '_niter_', num2str(n_iter), '_f_', num2str(f1), 'hz'];

for v = 1:length(values)
    value = values(v);

    [c_true, cmin, cmax] = GetTrueVelocity(nz, value, analysis_to_make);
   
    if(analysis_to_make == 2)
        plot(z, c_true, 'Color', colors{v}, 'LineWidth', linewidth, 'HandleVisibility','off'); 
    else
        if(v == 1)
            plot(z, c_true, ['-', 'b'], 'LineWidth', linewidth, 'HandleVisibility','on'); 
            lbl{1} = 'True velocity';
        end
    end
 
    source=SourceTerm(nt,f1,dt,print);
    p0 = zeros(nz,1);
    p1 = p0;
    p0(zs) = source(1);
    p1(zs) = source(2);
    
    [p_true, p_all] = ForwardProblem(nz, nt, source, c_true, zs, dx, dt, p0, p1, print);
    p_obs=p_true+noise_on_observed_pressure*randn(nt, 1);
   
    c = GetInitialVelocity(c_true, value, analysis_to_make);
    
    if(analysis_to_make == 3)
        plot(z,c, ':',' Color', colors{v}, 'LineWidth', linewidth, 'HandleVisibility','off')
    end

    fprintf('\n**********************\nRunning Algo Descent : peak frequency f= %d\n',f1);
 
    if(follow_optimisation)
        fprintf('it\t obj\t\t norm(d)\t step\n');
    end

    for i=1:n_iter 
        [objfun,g]= CostFunc_FWI(c);

        if(analysis_to_make==1)
            alpha = value/norm(g);
        else
            alpha = 100/norm(g);
        end

        d = -g;    
        c = c + alpha*d;    

        if(follow_optimisation == 1)
            fprintf('%d\t %.3e\t %.3e\t %.3e \n',i,objfun,norm(g), alpha);
            figure(3)
            plot(z,c_true,'b',z,c,'r');legend('True Velocity','Inverted Velocity');
            xlabel('Depth (m)');ylabel('Velocity (m/s)');
            title('Low frequency inversion');
            axis([0 1500 1000 3000]);
            drawnow;
        end

        if(norm(g) < tolg)
            fprintf("\n========== Stopping criteria met ==========\n")
            break;
        end
    end
    plot(z,c, 'Color', colors{v}, 'LineWidth', linewidth);
    if(analysis_to_make==0)
        lbl{v+1} = ['$c_0$=', num2str(value), ' m/s'];
    elseif(analysis_to_make==1)
        lbl{v+1} = ['$s_k \|g_k\| = $ ', num2str(value)];  
    elseif(analysis_to_make==2)
        lbl{v} = ['shift=', num2str(value * dx), ' m'];
    elseif(analysis_to_make==3)
        lbl{v} = ['span=', num2str(value)];
    end
    
end

axis([0 1500 1000 3000]);
grid minor
xl = xlabel('Depth (m)');
yl = ylabel('Velocity (m/s)');
set(yl, 'Interpreter', 'latex', 'fontsize', fontsize , 'LineWidth', linewidth);
set(xl, 'Interpreter', 'latex', 'fontsize', fontsize , 'LineWidth', linewidth);
leg = legend(lbl);
set(leg, 'Interpreter', 'latex', 'fontsize', legend_fontsize , 'LineWidth', linewidth, 'Location', 'southeast', 'NumColumns',2); % southeast bestoutside
clear print
print(4, titles, '-depsc')


% ---------------- TO PLOT PRESSURE ------------------- %
if(plot_pressure == 1)
    clear print
    eta_list = noise_on_observed_pressure;
    lbl = cell(1,length(eta_list) + 1);
    plot(t, p_true,'k', 'LineWidth', 1.2);
    lbl{1} = "True Pressure";
    for i=1:length(eta_list)
        eta = eta_list(i);
        p_obs = p_true+eta*randn(nt, 1);
        plot(t, p_obs, 'LineWidth', 1.2);
        lbl{i+1} = ['$\eta$=', num2str(eta)];
    end
    leg = legend(lbl);
    xl=xlabel('Depth (m)');
    yl=ylabel('Pressure (Pa)');
    set(yl, 'Interpreter', 'latex', 'fontsize', fontsize , 'LineWidth', linewidth);clc
    set(xl, 'Interpreter', 'latex', 'fontsize', fontsize , 'LineWidth', linewidth);
    set(leg, 'Interpreter', 'latex', 'fontsize', legend_fontsize , 'LineWidth', linewidth, 'Location', 'southeast', 'NumColumns', 1); % southeast bestoutside
    
    axis([0, nz*dt, -30, 30]);
    clear print
    titles = ['p_true_vs_obs_f:', num2str(f1)];
    print(4, titles, '-depsc')
end