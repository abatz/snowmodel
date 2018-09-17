function [AET,DEF,RUNOFF,SOILS]=simplebucketmodel(TMEAN,PPT,PET,AWC,SnowMelt,Sublimation,SnowDepth,tax,soil);
% PPT, PET, AWC should be in mm
% Sublimation is the water-flux already account for in the snow model
% PPT should be liquid precipitation plus snowmelt
% TMEAN in C
% tax is the fraction of total precip/snow melt that goes directly to runoff,
% nominally this is 0.05
% soil is initial conditions (can be generated by spinup), units of mm

if nargin==7  % if no initialization for soil moisture, start with 1/2 full
    soil=AWC/2;
end

AET=NaN*ones(size(PPT));DEF=AET;RUNOFF=AET;SOILS=AET;
% first let's correct PET for conditions when temperatures are cold, this 
MF=runsnow(TMEAN+273.15,1);
PET=PET.*(1-MF);
PPT=PPT.*(1-MF);
clear TMEAN

PET=PET-Sublimation;
% don't let PET<0,
f=find(PET<0);PET(f)=0;
PPT=PPT+SnowMelt;
INIT_RO=tax*PPT;
PPT=PPT*(1-tax);
lastsoil=soil;
for j=1:length(PPT)
    
deltasoil=PPT(j)-PET(j);

if deltasoil<0
    if deltasoil>lastsoil
        deltasoil=lastsoil;  % can't ask for more than is available in soil
    end
    if SnowDepth(j)<1 % this means we can take from soil column
% can not drain more than is in soil!!!
    drainsoil=-deltasoil.*(1-exp(-lastsoil./AWC));
    AET(j)=drainsoil+PPT(j);
    DEF(j)=PET(j)-AET(j);
    SOILS(j)=lastsoil-drainsoil;
    RUNOFF(j)=0+INIT_RO(j);
    else
        % this is where there is snowpack and deficit after sublimaton
        AET(j)=0; % could equal sublimation, but this is NOT used by plants
        DEF(j)=PET(j)-Sublimation(j);
        lastsoil=lastsoil; % soil is offlimits when snow is present!
        RUNOFF(j)=0+INIT_RO(j);
    end
else
    DEF(j)=0;
    AET(j)=PET(j);
    surplus=PPT(j)-AET(j);
    fillsoil=AWC-lastsoil;
    if fillsoil==0
        RUNOFF(j)=PPT(j)-AET(j)+INIT_RO(j);
        lastsoil=AWC;
        SOILS(j)=AWC;
    else
        if surplus>fillsoil
            lastsoil=AWC;
            SOILS(j)=AWC;
            RUNOFF(j)=surplus-fillsoil+INIT_RO(j);
        else
            lastsoil=lastsoil+surplus;
            SOILS(j)=lastsoil;
        end
    end
end
if isnan(AET(j))==1
    j
end

end
