function [lw]=longwave(emissivity,temp);

    SBconstant = 4.89E-06;%		#	[kJ m-2 K-4 d-1]
    tempK = temp + 273.15;%	#	[degrees K]
    lw=emissivity.*SBconstant.*tempK.^4;