function [FT,lMtilda,vMtilda,FMvtilda,FMpe] = RigidTendonForceVectorized2(a,lMT,vMT,params)

% Solve tendon force ODE for tendon force time history using numerical
% integration
% The equations in this function were generated by Autolev input file
% TendonForceODE.al

% Build up first order tendon force ODE dFT/dt = k*vT
% kT = b1*(FMo/lTs)-b2*((FMo*lTs)/(lMo^2))
% lM = sqrt((lMo*sin(alphao))^2+(lMT-lTs-FT/kT)^2)
% alpha = asin(lMo*sin(alphao)/lM)
% lMtilda = lM/lMo
% FMltilda = c11*exp(-c21*abs(lMtilda-c31)^c41)+c12*exp(-c22*abs(lMtilda-c32)^c42)+c13*exp(-c23*abs(lMtilda-c33)^c43)+c5
% FMce = (FT/cos(alpha))-FMo*(lMtilda-d1)^d2
% FMvtilda = FMce/(a*FMo*FMltilda)
% vMtilda = e1*sinh(e2*FMvtilda-e3)-e4+e5*(e2*FMvtilda-e3)+e6*(e2*FMvtilda-e3)^3
% vM = vMtilda*vMmax
% vT = vMT-vM/cos(alpha)
% dFT/dt = kT*vT

% Extract muscle-tendon model parameter values
FMo = ones(size(a,1),1)*params(1,:);
lMo = ones(size(a,1),1)*params(2,:);
lTs = ones(size(a,1),1)*params(3,:);
alphao = pi/180.0.*ones(size(a,1),1)*params(4,:); % Converts to radians
vMmax = ones(size(a,1),1)*params(5,:);
% b1 = ones(size(a,1),1)*params(6,:);
% b2 = ones(size(a,1),1)*params(7,:);
c11 = ones(size(a,1),1)*params(8,:);
c21 = ones(size(a,1),1)*params(9,:);
c31 = ones(size(a,1),1)*params(10,:);
c41 = ones(size(a,1),1)*params(11,:);
c12 = ones(size(a,1),1)*params(12,:);
c22 = ones(size(a,1),1)*params(13,:);
c32 = ones(size(a,1),1)*params(14,:);
c42 = ones(size(a,1),1)*params(15,:);
c13 = ones(size(a,1),1)*params(16,:);
c23 = ones(size(a,1),1)*params(17,:);
c33 = ones(size(a,1),1)*params(18,:);
c43 = ones(size(a,1),1)*params(19,:);
c5 = ones(size(a,1),1)*params(20,:);
d1 = ones(size(a,1),1)*params(21,:);
d2 = ones(size(a,1),1)*params(22,:);
% Use for rigid tendon model (Fm = f(vm))
e1 = ones(size(a,1),1)*params(29,:);
e2 = ones(size(a,1),1)*params(30,:);
e3 = ones(size(a,1),1)*params(31,:);
e4 = ones(size(a,1),1)*params(32,:);
e5 = ones(size(a,1),1)*params(33,:);
e6 = ones(size(a,1),1)*params(34,:);
% Parameters for lM equation
f1 = ones(size(a,1),1)*params(35,:);
f2 = ones(size(a,1),1)*params(36,:);
% Parameters for alpha curve
g1 = ones(size(a,1),1)*params(37,:);
g2 = ones(size(a,1),1)*params(38,:);
% Parameters for FT curve (not used in rigid tendon model)
% h1 = ones(size(a,1),1)*params(39,:);
% h2 = ones(size(a,1),1)*params(40,:);

% keyboard
% Calculate tendon forces for a rigid tendon model
lMTtilda = lMT./lTs;
lM = sqrt((lMo.*sin(alphao)).^2+(lTs.*(0.5.*(f1.*log(cosh((lMTtilda-1)./f1))+lMTtilda)-f2)).^2);
lMtilda = lM./lMo;
arg = lMo.*(sin(alphao))./lM;
alpha = asin(-0.5.*(g1.*log(cosh((arg-1)./g1))-arg)+g2);
FMltilda = c11.*exp(-c21.*((lMtilda-c31).^2).^(c41/2))+c12.*exp(-c22.*((lMtilda-c32).^2).^(c42/2))+c13.*exp(-c23.*((lMtilda-c33).^2).^(c43/2))+c5;
vM = vMT.*cos(alpha);
vMtilda = vM./vMmax;
FMvtilda = -e1.*asinh(-e2.*vMtilda-e3)+e4-e5.*(-e2.*vMtilda-e3)+e6.*(-e2.*vMtilda-e3).^3;
FMce = a.*FMo.*FMltilda.*FMvtilda;
FMpe = FMo.*((lMtilda-d1).^2).^(d2/2);
FM = FMce+FMpe;
FT = FM.*cos(alpha);

return