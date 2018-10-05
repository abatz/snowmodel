function [Correction]=StabilityCorrection(Z,d,Tsurf,Tair,Wind);

RiCr = 0.2;           %/* Critical Richardson's Number */
Correction = 1.0;
CONST_TKFRZ=273.15;
CONST_G=9.80616;


%Calculate the effect of the atmospheric stability using a Richardson Number approach */

if Tsurf ~= Tair %       /* Non-neutral conditions */
    Ri = CONST_G * (Tair - Tsurf) * (Z ) /(((Tair + CONST_TKFRZ) + (Tsurf + CONST_TKFRZ)) / 2.0 * Wind^2);
    Ri=min(Ri,RiCr);
    if Ri>0
      Correction = (1 - Ri / RiCr)^2; 
    else
      Correction = real(sqrt(1 - 16 * Ri));
    end
end