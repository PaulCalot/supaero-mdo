%
%          Y. Diouane
%          BE: Problemes Inverses
%          ISAE-SUPAERO
%
% Question 1:
%
clear all; close all;
global nt dt dx zs source p_obs cmin cmax
rand('seed',1)
randn('seed',1)
nt = 1500;
%To avoid numerical dispersion and instability, dx < cmin/(5*fpeak) and dt < 0.6*dx/cmax.
figure(3)
hold on
values = [2 3 5 7 8];
lbl = cell(size(values));
unit = 'm';
quant = 'z_s';
quant_lab = 'z_source';

for v = 1:length(values)
    value = values(v);
    
    dx = 5;
    dt = 0.001;
    nz = 301;
    c_true = 2000*ones(nz,1);
    c_true(51:80) = 1600;
    c_true(121:150) = 2500;
    c_true(200:210) = 1800;
    %c_true = c_true + 200*value*rand(size(c_true));
    
    cmin=min(c_true);
    cmax=max(c_true);
    zs = value; % The source position

    f1 = 5; % The peak frequency 
    print=0;
    source=SourceTerm(nt,f1,dt,print);
    p0 = zeros(nz,1);
    p1 = p0;
    p0(zs) = source(1);
    p1(zs) = source(2);
    [p_true,p_all] = ForwardProblem(nz,nt,source,c_true,zs,dx,dt,p0,p1,print);
    p_obs = p_true ; % possibly with noise
    %
    n = 100;
    t = linspace(-2,2,n);  % t \in [-2, 2]
    c0 = c_true;
    c1 = c0 + 0.5*c0;
    %
    % define vect_f such as vect_f(i) = obj_fun(c+ t(i) dc)
    %

    for i =1:n
        
        vect_f(i) = CostFunc_FWI(c0 + t(i)*(c1-c0)); %#ok<*SAGROW>
    end

plot(t, vect_f,'LineWidth',2);
lbl{v} = ['$',quant, ' = ',num2str(value),' ',unit,'$'];
end

xl = xlabel({'t'});
set(xl, 'Interpreter', 'latex', 'fontsize', 16 , 'LineWidth', 1.5);
yl = ylabel({'Obj. Func.'});
set(yl, 'Interpreter', 'latex', 'fontsize', 16 , 'LineWidth', 1.5);
leg = legend(lbl);
set(leg, 'Interpreter', 'latex', 'fontsize', 14, 'Location', 'bestoutside')
ax = gca;
ax.FontSize = 16;
grid minor
titles = ['Quest1 for var ', quant_lab];
clear print
print(3,titles,'-depsc');