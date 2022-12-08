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
%
dx = 5;
dt = 0.001;
nz = 301;
c_true = 2000*ones(nz,1);
c_true(51:80) = 1600;
c_true(121:150) = 2500;
c_true(200:210) = 1800;
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
[p_true,p_all] = ForwardProblem(nz,nt,source,c_true,zs,dx,dt,p0,p1,print);
p_obs=p_true; % possiblly with noise
%
n = 100;
t = linspace(-2,2,n);  % t \in [-2, 2]
%
% define vect_f such as vect_f(i) = obj_fun(c+ t(i) dc)
%
for i =1:n
    %
    % To Do
    %
    vect_f(i)=
end
figure(3);
plot(t, vect_f,'k','LineWidth',2);
xlabel({'t'});
ylabel({'Obj. Func.'});