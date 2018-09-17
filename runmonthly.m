% load in sample file with requisite monthly climate data
% contains ppt (mm), srad (W/m2), tmax(K), tmin(K), sph(kg/kg), vs (m/s)

load sampleclimateinputs
elevation=2000; % you can change this (m)
SnowMelt=NaN*ones(size(ppt));SnowWaterEq=SnowMelt;
maxsolar=potential_solar(lat,[1:365],elevation); 
for i=1:12
    d1=datenum(2013,i,1)-datenum(2013,1,0);
    d2=datenum(2013,i+1,1)-1-datenum(2013,1,0);
    maxsun(i)=nanmean(maxsolar(d1:d2));
end

% assume that minimum solar insolation is 35% of maximum insolation for
% clear sky conditions
minsun=maxsun*.35;

cloudiness=1-(srad-repmat(minsun',[1 40]))./(repmat(maxsun'-minsun',[1 40]));
srad=srad*.0864*1000;

tmax=tmax-273.15;
tmin=tmin-273.15;
tdew=dewpoint(sph,elevation);
tdew=tdew-273.15;
[SnowMelt,SnowWaterEq]=runsnowmelt_monthly(cloudiness(:),tmax(:),tmin(:),ppt(:),srad(:),tdew(:),vs(:),elevation);
[SnowMelt1,SnowWaterEq1,X,Sublimation]=runsnowmelt_monthly_wlatent(cloudiness(:),tmax(:),tmin(:),ppt(:),srad(:),tdew(:),vs(:),elevation);
figure(1);
subplot(2,1,1)
plot(1979:1/12:2018.99,SnowWaterEq(:),'k');hold on;plot(1979:1/12:2018.99,SnowWaterEq1(:),'b');legend({'SWE';'SWE_{sublimation}'})
ylabel('SWE (m)');axis tight
subplot(2,1,2)
bar(12*[nanmean(SnowMelt1(:)) nanmean(max(0,Sublimation(:)))]);set(gca,'XTickLabel',{'Snow Melt';'Sublimation'});
ylabel('Annual Mean Loss of SWE (m)');
