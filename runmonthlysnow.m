lat=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/rmax/rmax_2017.nc','lat');
lon=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/rmax/rmax_2017.nc','lon');
y=70;x=190;
sph=matfile('/Users/drthunder/metdata_sph_CONUS_monthly.mat');
tmin=matfile('/Users/drthunder/metdata_tmin_CONUS_monthly.mat');
tmax=matfile('/Users/drthunder/metdata_tmax_CONUS_monthly.mat');
ppt=matfile('/Users/drthunder/metdata_ppt_CONUS_monthly.mat');
vs=matfile('/Users/drthunder/metdata_vs_CONUS_monthly.mat');
srad=matfile('/Users/drthunder/metdata_srad_CONUS_monthly.mat');
SnowMelt=NaN*ones(388,400,468);SnowWaterEq=SnowMelt;
for y=13:400;
maxsun=max(max(squeeze(srad.data(y,:,:,:)),[],1),[],3);
minsun=min(min(squeeze(srad.data(y,:,:,1:38)),[],1),[],3);
    for x=1:400
pptdata=squeeze(ppt.data(y,x,:,:));
if ~isnan(pptdata(3))
tmaxdata=squeeze(tmax.data(y,x,:,:));
tmindata=squeeze(tmin.data(y,x,:,:));
sphdata=squeeze(sph.data(y,x,:,:));
vsdata=squeeze(vs.data(y,x,:,:));
solardata=squeeze(srad.data(y,x,:,:));

cloudiness=1-(solardata-repmat(minsun',[1 39]))./(repmat(maxsun'-minsun',[1 39]));
solardata=solardata*.0864*1000;

tmaxdata=tmaxdata-273.15;
tmindata=tmindata-273.15;
tdew=dewpoint(sphdata,1000);
tdew=tdew-273.15;
tmaxdata(468)=tmaxdata(467);tmindata(468)=tmindata(467);cloudiness(468)=cloudiness(467);
vsdata(468)=vsdata(467);tdew(468)=tdew(467);solardata(468)=solardata(467);
[SnowMelt(y-12,x,:),SnowWaterEq(y-12,x,:)]=runsnowmelt_monthly(cloudiness(:),tmaxdata(:),tmindata(:),pptdata(:),solardata(:),tdew(:),vsdata(:),1000);
end
    end
    y
end
