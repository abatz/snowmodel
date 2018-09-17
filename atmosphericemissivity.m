function [AE]=atmosphericemissivity(temp,cloud,tdmean,elev);
%AtmosphericEmissivity <- function (airtemp, cloudiness, vp = NULL, opt = "linear") 
%	# The emissivity of the atmsophere [-]
%		# airtemp:	air temperature [C]
%		# cloudiness : fraction of the sky covered in clouds [-]
%		# vp : Vapor Pressure :  [kPa]
vp=calcVP(tdmean,elev);
vp=vp/10;

AE= 1.24 * (((vp * 10)./(temp + 273.15)).^(1/7) .* (1 - 0.84 * cloud) + 0.84 * cloud);
