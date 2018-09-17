function [vp]=calcVPD(tdew,Z);
% SATVAP: computes saturation vapor pressure
% q=satvap(Ta) computes the vapor pressure at satuation at air
% temperature Ta (deg C). From Gill (1982), Atmos-Ocean Dynamics, 606.
%
%    INPUT:   Ta- air temperature  [C]
%             Z - elevation
%
%    OUTPUT:  vp  - saturation vapour pressure  [mb]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
% 8/27/98: version 1.1 (corrected by RP)
% 8/5/99: version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% send data in degrees C

P=1013.25*(1-.0001*Z);

% solve for vapor pressure at maximum temp

ew=power(10,((0.7859+0.03477*tdew)./(1+0.00412*tdew)));
fw=1 + 1e-6*P.*(4.5+0.0006*tdew.^2);
ew_tdew=fw.*ew;

vp=ew_tdew;

