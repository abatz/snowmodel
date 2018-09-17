% need to load solar, tmax, tmin, tdew, vs, solar, ppt, elevation
srad=load800(1,1979);
vs=load800(2,1979);
tmax=load800(3,1979);
tmin=load800(4,1979);
tdmean=load800(6,1979);
[ppt,lon,lat,elev]=load800(5,1979);
cloudiness=load800(7,1979);
srad=srad*.0864*1000;
SnowMelt=single(NaN*ones(size(tmax)));
SnowWaterEq=SnowMelt;
StartSWE=zeros(3105,7025);
for y=1:3105;
    parfor x=1:7025
[SnowMelt(y,x,:),SnowWaterEq(y,x,:),SnowDepth(y,x)]=runsnowmelt_monthly(squeeze(cloudiness(y,x,:)),squeeze(tmax(y,x,:)),squeeze(tmin(y,x,:)),squeeze(ppt(y,x,:)),squeeze(srad(y,x,:)),squeeze(tdmean(y,x,:)),squeeze(vs(y,x,:)),elev(y,x));
end
    end
save(['SWE_',num2str(1978+i)],'-v7.3','SnowMelt','SnowWaterEq','SnowDepth');
StartSWE=SnowWaterEq(:,:,12);
clear SnowMelt SnowWaterEq cloudiness tmax tmin ppt srad tdmean vs 

for i=1:38
srad=load800(1,1979+i);
vs=load800(2,1979+i);
tmax=load800(3,1979+i);
tmin=load800(4,1979+i);
tdmean=load800(6,1979+i);
[ppt,lon,lat,elev]=load800(5,1979+i);
cloudiness=load800(7,1979+i);
srad=srad*.0864*1000;
SnowMelt=single(NaN*ones(size(tmax)));
SnowWaterEq=SnowMelt;
for y=1:3105;
    parfor x=1:7025
[SnowMelt(y,x,:),SnowWaterEq(y,x,:),SnowDepth(y,x)]=runsnowmelt_monthly(squeeze(cloudiness(y,x,:)),squeeze(tmax(y,x,:)),squeeze(tmin(y,x,:)),squeeze(ppt(y,x,:)),squeeze(srad(y,x,:)),squeeze(tdmean(y,x,:)),squeeze(vs(y,x,:)),elev(y,x),StartSWE(y,x),SnowDepth(y,x));
end
    end
save(['SWE_',num2str(1979+i)],'-v7.3','SnowMelt','SnowWaterEq','SnowDepth');
StartSWE=SnowWaterEq(:,:,12);
clear SnowMelt SnowWaterEq cloudiness tmax tmin ppt srad tdmean vs 
i
end
