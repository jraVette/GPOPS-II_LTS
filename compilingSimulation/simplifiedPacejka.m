function [Fx, Fy, muX, muY, FxMax, FyMax, kappa_n, alpha_n, rho, efficiency] = simplifiedPacejka(Fz,alpha,kappa,coeff)
%This function is a simplifed pacejka model explained in several places in
%the literature.  Appendix 1 of:
%
%G. Perantoni and D. J. Limebeer, ?Optimal control for a formula one car
%with variable parameters,? Vehicle System Dynamics, vol. 52, no. 5, pp.
%653?678, 2014.
%
%Is where this function came from
%
%INPUTS:
%    fz - tire vertical load in N                              class double
%    alpha - tire slip angle in rad                            class double
%    kappa - slip ratio
%    coeff - coeff structure of the tire model 

%Creation 4 Aug 2015 - Jeff Anderson

    d2r = pi/180;


%If not tire is passed in, then use nominal timre model from paper
% if nargin<4
%     coeff.Fz1     = 2000;     %[N] Reference Load 1
%     coeff.Fz2     = 6000;     %[N] Reference Load 2
%     coeff.muX1    = 1.75;     %Peak longitudinal friction coeff at load 1
%     coeff.muX2    = 1.40;     %Peak longitidnal friction coeff at load 2
%     coeff.kappa1  = 0.11;     %slip coeff for friction peak at load 1
%     coeff.kappa2  = 0.10;     %Slip coeff for friction peak at load 2
%     coeff.muY1    = 1.80;     %Peak lat friction coeff at load 1
%     coeff.muY2    = 1.45;     %Peak lat friction coeff at load 2
%     coeff.alpha1  = 9*d2r;    %Slip angle for friction peak at load 1
%     coeff.alpha2  = 8*d2r;    %Slip angle for friction peak at load 1
%     coeff.Qx      = 1.9;      %Longutindal shape factor
%     coeff.Qy      = 1.9;      %Lateral shape factor
% end

%Convet to rad
% alpha = alpha*d2r;

% %Peak mu's
muX_max  = interp1([coeff.Fz1 coeff.Fz2],[coeff.muX1   coeff.muX2],  Fz,'linear','extrap');
muY_max  = interp1([coeff.Fz1 coeff.Fz2],[coeff.muY1   coeff.muY2],  Fz,'linear','extrap');
kappaMax = interp1([coeff.Fz1 coeff.Fz2],[coeff.kappa1 coeff.kappa2],Fz,'linear','extrap');
alphaMax = interp1([coeff.Fz1 coeff.Fz2],[coeff.alpha1 coeff.alpha2],Fz,'linear','extrap');

% muX_max = (Fz - coeff.Fz1).*(coeff.muX2 - coeff.muX1)./(coeff.Fz2 - coeff.Fz1) + coeff.muX1;
% muY_max = (Fz - coeff.Fz1).*(coeff.muY2 - coeff.muY1)./(coeff.Fz2 - coeff.Fz1) + coeff.muY1;
% kappaMax = (Fz - coeff.Fz1).*(coeff.kappa2 - coeff.kappa1)./(coeff.Fz2 - coeff.Fz1) + coeff.kappa1;
% alphaMax = (Fz - coeff.Fz1).*(coeff.alpha2 - coeff.alpha1)./(coeff.Fz2 - coeff.Fz1) + coeff.alpha1;

%Peak forces possible
FxMax = muX_max.*Fz;
FyMax = muY_max.*Fz;


%Normalized slip
kappa_n = kappa./kappaMax;
alpha_n = alpha./alphaMax;

rho = sqrt(alpha_n.^2 + kappa_n.^2);

%Friction coeff
Sx = pi/(2*atan(coeff.Qx));
Sy = pi/(2*atan(coeff.Qy));

muX = muX_max.*sin(coeff.Qx*atan(Sx*rho));
muY = muY_max.*sin(coeff.Qy*atan(Sy*rho));


Fx = muX.*Fz.*kappa_n./(rho+1e-16);
Fy = muY.*Fz.*alpha_n./(rho+1e-16); %use eps to eliemitae nan is /0


efficiency = sqrt((Fx./FxMax).^2 + (Fy./FyMax).^2);

    