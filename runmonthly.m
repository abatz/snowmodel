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

% run normal snowmelt model w/o contribution from sublimation
[SnowMelt,SnowWaterEq,X,Sublimation,SnowTemp,E,Energy,Albedo]=runsnowmelt_monthly_sublimation(cloudiness(:),tmax(:),tmin(:),ppt(:),srad(:),tdew(:),vs(:),elevation,0);
% run normal snowmelt model w/contribution from sublimation
[SnowMelt1,SnowWaterEq1,X1,Sublimation,SnowTemp1,E1,Energy1,Albedo]=runsnowmelt_monthly_sublimation(cloudiness(:),tmax(:),tmin(:),ppt(:),srad(:),tdew(:),vs(:),elevation,1);

figure(1);clf;
subplot(2,1,1)
plot(1979:1/12:2018.99,SnowWaterEq(:),'k');hold on;plot(1979:1/12:2018.99,SnowWaterEq1(:),'b');legend({'SWE';'SWE_{sublimation}'})
ylabel('SWE (m)');axis tight
subplot(2,1,2)
bar(12*[nanmean(SnowMelt1(:)) nanmean(max(0,Sublimation(:)))]);set(gca,'XTickLabel',{'Snow Melt';'Sublimation'});
ylabel('Annual Mean Loss of SWE (m)');

% now run a simple water balance model
figure(2);clf;

f=find(SnowWaterEq1<.001);Albedo(f)=0.23;
srad2=srad/.0864/1000;
[PET]=monthlyPET(srad2, tmax,tmin, vs,lat,elevation,Albedo,tdew,40);

% available water holding content (AWC) should come from observations
 AWC=150;  
 TMEAN=tmax/2+tmin/2;
 tax=0.05;
 soil=100; % initialize soil moisture
 TMEAN=TMEAN(:);ppt=ppt(:);PET=PET(:);
 [AET,DEF,RUNOFF,SOILS]=simplebucketmodel(TMEAN,ppt,PET,AWC,SnowMelt1*1e3,Sublimation*1e3,SnowWaterEq1*1e3,tax,soil);
 AET=reshape(AET,size(tmax));
 DEF=reshape(DEF,size(tmax));
 SOILS=reshape(SOILS,size(tmax));
 RUNOFF=reshape(RUNOFF,size(tmax));
 subplot(2,2,1)
 plot(nanmean(SOILS,2));
 hold on
 plot(nanmean(AET,2),'r');
 plot(nanmean(DEF,2),'k');
 legend({'SOIL';'AET';'DEF'})
 ylabel('mm');
 
  subplot(2,2,2)
 plot(nanmean(reshape(ppt+SnowMelt1*1e3,12,40),2));
 hold on
 plot(nanmean(RUNOFF,2),'r');
 plot(nanmean(reshape(ppt,12,40),2),'k');
 legend({'Rain+Melt';'Runoff';'PPT'})
  ylabel('mm');

 
 subplot(2,1,2);
 plot(1979:2017,zscore(nanmean(SOILS(6:9,1:39),1)));
 hold on;
 plot(1979:2017,zscore(nansum(AET(:,1:39),1)));
 hold on;
 plot(1979:2017,zscore(nansum(DEF(:,1:39),1)));
 legend({'SOIL';'AET';'DEF'})
 ylabel('\sigma');