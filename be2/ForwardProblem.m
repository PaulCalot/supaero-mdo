function [p_calc,p_all] = ForwardProblem(nz,nt,source,c,zs,dx,dt,p0,p1,print)
%
%          Y. Diouane
%          BE: Problemes Inverses
%          ISAE-SUPAERO
%  
%To avoid numerical dispersion and instability, dx < cmin/(5*fpeak) and dt < 0.6*dx/cmax.
p = zeros(size(p0));
p_all = zeros(nz,nt);
p_calc = zeros(nt,1);
p_calc(1) = p0(zs);
p_calc(2) = p1(zs);
p_all(:,1) = p0;
p_all(:,2) = p1;
for it=3:nt
    for iz=2:nz-1
        % CFL Condition : keep alpha small
        alpha = c(iz)*dt/dx;
        alpha = alpha*alpha;
        p(iz) = 2*(1-alpha)*p1(iz) + alpha*(p1(iz+1)+p1(iz-1)) - p0(iz);
    end
    p(zs) = p(zs) + source(it);
    p_calc(it) = p(zs);
    p0 = p1;
    p1 = p;
    p_all(:,it) = p;
end
 if print == 1
    step_t=4;
    for it=1:step_t:nt
    plot(p_all(:,it));legend('The wave form');
    xlabel('Depth (m)');ylabel('Wave form');
    axis([0 nz -30 30]);
    drawnow
    end
end
