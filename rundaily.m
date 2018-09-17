% load data directly from NKN thredds server for a given lat/lon
% this runs the daily snowmodel
lat=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/rmax/rmax_2017.nc','lat');
lon=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/rmax/rmax_2017.nc','lon');
sph=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/sph/sph_2017.nc','specific_humidity',[63 190 1],[1 1 365],[1 1 1]);
tmax=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/tmmx/tmmx_2017.nc','air_temperature',[63 190 1],[1 1 365],[1 1 1]);
tmin=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/tmmn/tmmn_2017.nc','air_temperature',[63 190 1],[1 1 365],[1 1 1]);
ppt=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/pr/pr_2017.nc','precipitation_amount',[63 190 1],[1 1 365],[1 1 1]);
vs=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/vs/vs_2017.nc','wind_speed',[63 190 1],[1 1 365],[1 1 1]);
solar=ncread('http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/srad/srad_2017.nc','surface_downwelling_shortwave_flux_in_air',[63 190 1],[1 1 365],[1 1 1]);
ppt=squeeze(ppt);
tmax=squeeze(tmax);
tmin=squeeze(tmin);
sph=squeeze(sph);
vs=squeeze(vs);
solar=squeeze(solar)*.0864*1000;
cloudiness=ones(size(solar));

tmax=tmax-273.15;
tmin=tmin-273.15;
tdew=dewpoint(sph,1000);
tdew=tdew-273.15;
[SnowMelt,SnowWaterEq]=runsnowmelt(cloudiness,tmax,tmin,ppt,solar,tdew,vs,elev);
