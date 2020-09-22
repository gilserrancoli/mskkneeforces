function [b,c,d,ev,eF,f,g,h] = GetNormalizedParameterValues(problem)

% % Tendon model parameter values for tendon stiffness
% kT = (FMo/lTs)*(b1-b2*(lTs/lMo)^2) = b1*(FMo/lTs)-b2*((FMo*lTs)/(lMo^2))
% where kT is non-normalized tendon stiffness
% Note that making b2 = 0.4 will cause the plantar flexor muscles to have
% a more compliant tendon than when b2 = 0, consistent with plantar flexor
% tendon stiffness values reported in the literature, while having little
% influence on the tendon stiffness values for other muscles
% Using a b1 value of 37.5 is consistent with the Zajac data
% Using a b1 value of 20 is consistnet with the Millard paper (and
% experimental data) Note that the choice of the b1 value affects the
% choice of the h2 parameter for the FT equation for the compliant tendon

b1 = 20; 
b2 = 0;
b = [b1; b2];

% Muscle model paramter values for normalized active force-length curve
% FMactive = sum(c1i*exp(-c2i*abs(lM-c3i).^c4i))+c5 for i = 1 to 3
% FMactive is normalized muscle force and lM is normalized muscle length
c11 = 9.404711698062630e-01;
c21 = 1.073124043936821e+01;
c31 = 9.871922298750729e-01;
c41 = 1.967912928195565e+00;
c12 = 3.553707328939712e-01;
c22 = 4.446159579601148e+01;
c32 = 6.499995398513716e-01;
c42 = 2.061727424469695e+00;
c13 = 3.319108437630654e-01;
c23 = 2.990238558320256e+01;
c33 = 1.432634023963042e+00;
c43 = 2.524739871950831e+00;
c5 = 5.000000000000000e-02;
c = [c11; c21; c31; c41; c12; c22; c32; c42; c13; c23; c33; c43; c5];

% Muscle model parameter values for normalized passive force-length curve
% FMpassive = abs(lM-d1).^d2;
% where FMpassive is normalized muscle force and lM is normalized muscle length
d1 = 0.6;
d2 = 5;
d = [d1; d2];

% Muscle model parameter values for normalized force-velocity curves
% Normalized muscle velocity vM as a function normalized muscle force FM
% vM = e1*sinh(x)-e4+e5*x+e6*x^3 where x = e2*FM-e3
% where vM is normalized muscle velocity and FM is normalized muscle force
switch problem
    
    case {'max','submax'}
        % vM = f(FM) fitted using muscle force data from maximum activation
        % biological benchmark problem
        e1 = 1.810014757034373e-05;
        e2 = 1.228831369367371e+01;
        e3 = 1.104066589257847e+01;
        e4 = 1.312235742558771e-02;
        e5 = 3.500979685877418e-03;
        e6 = 2.854136006923058e-04;
        ev = [e1; e2; e3; e4; e5; e6];
        
    otherwise
        % vM = f(FM) fitted using default OpenSim parameter values reported in
        % https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1ForceVelocityInverseCurve.html
        e1 = 5.050895036092566e-06;
        e2 = 1.385120877069190e+01;
        e3 = 1.246603410122264e+01;
        e4 = 4.792117529500456e-03;
        e5 = 6.529188684572741e-03;
        e6 = 1.295085984254243e-04;
        ev = [e1; e2; e3; e4; e5; e6];

end
        
% Normalized muscle force FM as a function of normalized muscle velocity vM
% FM = -e1*asinh(x)+e4-e5*x+e6*x^3 where x = -e2*vM-e3
% (or equivalently FM = -e1*log(x+sqrt(x^2+1))+e4-e5*x+e6*x^3)
% where FM is normalized muscle force and vM is normalized muscle velocity
% Note that FM = f(vM) curves produced by the parameter values below are
% identical to the inverted vM = f(FM) curves produced by the parameter
% values above
switch problem
    
    case {'max','submax'}
        % FM = f(vM) fitted using muscle force data from maximum activation
        % biological benchmark problem
        e1 = 1.600865120640679e-01;
        e2 = 1.157830112768191e+02;
        e3 = 1.517542572900223e+00;
        e4 = 8.982647358657909e-01;
        e5 = 1.061737623184449e-03;
        e6 = 6.663738886311849e-08;
        eF = [e1; e2; e3; e4; e5; e6];

    otherwise
        % Slope of 3
        % FM = f(vM) fitted using default OpenSim parameter values reported in
        % https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1ForceVelocityInverseCurve.html
        e1 = 0.708716646303580;
        e2 = 5.65697837083118;
        e3 = 0.191884073825674;
        e4 = 0.896633651436123;
        e5 = -0.165769379169328;
        e6 = -0.000623692878996660;
        eF = [e1; e2; e3; e4; e5; e6];
%         % Slope of 10
%         % FM = f(vM) fitted using default OpenSim parameter values reported in
%         % https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1ForceVelocityInverseCurve.html
%         e1 = 2.549488812131978e-01;
%         e2 = 4.371028583142105e+01;
%         e3 = 4.178636024718488e-01;
%         e4 = 8.990442378726738e-01;
%         e5 = -6.458505959106526e-03;
%         e6 = -5.084516555772948e-07;
%         eF = [e1; e2; e3; e4; e5; e6];

end

% Muscle model parameter values for muscle length equation in the rigid 
% tendon model that replaces (lMT-lTs) in the original equation: 
% lM = sqrt((lMo*sin(alphao))^2+(lMT-lTs)^2)
% with an analytic function that goes to zero when lMT < lTs.
% The new lM equation with the analytic function is:
% lM = sqrt((lMo*sin(alphao))^2+(lTs*(0.5*(f1*log(cosh((lMTtilda-1)/f1))+lMTtilda)-f2))^2)
% where lM is muscle length, lMo is optimal fiber length, lTs is tendon
% slack length, and lMTtilda = lMT/lTs, where lMT is muscle-tendon length
f1 = 0.05;
f2 = 0.482671320537530;
f = [f1; f2];

% Muscle model parameter values for pennation angle equation in the rigid 
% and compliant tendon model to prevent an issue where the input to the 
% asin function goes slightly higher than 1. The new equation is: 
% alpha = asin(-0.5*(g1*log(cosh((arg-1)/g1))-arg)+g2)
% where alpha is the pennation angle and arg is lMo*(sin(alphao))/lM
g1 = 0.005;
g2 = 0.498267132048600;
g = [g1; g2];

% Muscle model parameter values for tendon force equation in the compliant 
% tendon model that replaces (lT-lTs) in the original equation: 
% FT = kT.*(lT-lTs)
% with an analytic function that goes to zero when lT < lTs.
% The new FT equation with the analytic function is:
% FT = FMo*(0.5*kT*(h1*log(cosh((-1+lTtilda)/h1))+lTtilda)-h2)
% where FMo is maximum isometric strength, kT is a constant (either 37.5 or
% 20, see description of b parameters above), and lTilda is lT/lTs, where
% lT is tendon length and lTs is tendon slack length.
% Note that the value of h2 depends on the value of b1.
h1 = 0.000075;
if b1 == 37.5
    h2 = 18.7490252617773;
elseif b1 == 20
    h2 = 9.99948013961458;
else
    keyboard
    printf('h2 parameter needs to be calculated based on the kT value.\n');
end
h = [h1; h2];

% Note that when the parameter vector params is formed, it will contain
% muscle and tendon model parameter values for each muscle in the
% following order:
% params(1,:) = FMo;
% params(2,:) = lMo;
% params(3,:) = lTs;
% params(4,:) = alphao; % (in degrees)
% params(5,:) = vMmax;
% params = [params; b; c; d; ev; ef; f; g; h];

return