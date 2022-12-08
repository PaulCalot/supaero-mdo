function [varargout]=CostFunc_FWI(c,varargin)

% Evaluate objective function (and the gradient)
% Usage:       f = CostFunc_FWI(c)           evaluates function value only
%              g = CostFunc_FWI(c,'gradient')     evaluates gradient only
%          [f,g] = CostFunc_FWI(c)  evaluates function value and gradient
%
%          Y. Diouane
%          BE: Problemes Inverses
%          ISAE-SUPAERO
%  
global nt dx dt zs source p_obs cmin cmax
%%To avoid numerical dispersion and instability, 
% we impose that cmin <= c <= cmax
c=max(cmin,c);
c= min(cmax,c); 
nz=length(c);
%
% MODELING (FORWARD PROBLEM)
% Compute synthetic data (data) using the trial velocity c
%
p0 = zeros(nz,1);
p1 = p0;
zs = 5;
p0(zs) = source(1);
p1(zs) = source(2);
[p_calc,p_all]=ForwardProblem(nz,nt,source,c,zs,dx,dt,p0,p1,0);
fun = (1/2)*norm(p_calc-p_obs)^2;
%
% Back propagate the residual
%
dp = p_obs - p_calc;
dp = flipud(dp);
p0(zs) = dp(1);
p1(zs) = dp(2);
[ddata,dp_all] = ForwardProblem(nz,nt,dp,c,zs,dx,dt,p0,p1,0);
%
% Compute the gradient vector
%
grad = sum(p_all.*fliplr(dp_all),2)./(c.^3);
if nargin == 1
    % Compute objective function value
    if nargout== 2
        % Gradient is requested
        varargout{1} =  fun;
        varargout{2} =  grad;
    else
        varargout{1} =  fun;
    end
elseif(nargin==2)
    if(strcmp(varargin{1},'gradient'))
        % Only gradient is requested
        varargout{1}=grad;
    end
end