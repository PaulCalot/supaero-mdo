function [source]=SourceTerm(nt,fr,dt,print)
%
%          Y. Diouane
%          BE: Problemes Inverses
%          ISAE-SUPAERO
%  
npt=nt*dt;
t=(-npt/2):dt:npt/2;
rick1=(1-t .*t * fr^2 *pi^2  ) .*exp(- t.^2 * pi^2 * fr^2 ) ;
rick=rick1(round(nt/2)-round(1/fr/dt)+1:nt);
source = [rick'; zeros(nt-length(rick),1)];
if(print==1)
l=length(source);
figure(4)
plot([0:l-1]*dt,source);xlabel('Time (s)');
title([num2str(fr),' Hz Ricker Wavelet']);pause(1);
end
