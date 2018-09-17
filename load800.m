function [data,lon,lat,elev]=load800(var,year);

m=matfile('/home/abatz/sradratio800.mat');
sradratio=m.sradratio;
el2=matfile('/data/obs/obs/gridded/metdata/dem/metdata_elev');el2=el2.elev;
load /data/obs/obs/gridded/metdata/alt_gridmetlonlat
[lona,lata]=meshgrid(lon,lat);

load /data/obs/obs/gridded/prism/dem/dem_800mPRISM_US.mat
[lon,lat]=meshgrid(lon,lat);
switch var,

case 1,
sradm=matfile('/data/obs/obs/gridded/metdata/metdata_srad_CONUS_monthly.mat');
srad=sradm.data(:,:,:,year-1978);
for j=1:12
data(:,:,j)=interp2(lona,lata,srad(:,:,j),lon,lat);
end
clear srad vs
data=data.*sradratio;

case 2,
vsm=matfile('/data/obs/obs/gridded/metdata/metdata_vs_CONUS_monthly.mat');
vs=vsm.data(:,:,:,year-1978);
for j=1:12
data(:,:,j)=interp2(lona,lata,vs(:,:,j),lon,lat);
end


case 3,
if year<2010
tmax=matfile(['/data/obs/obs/gridded/prism/lt71m/800m_grids/tmax/tmax_800mPRISM_',num2str(year),'_CONUS.mat']);
data=tmax.data;
else
tmaxm=matfile('/data/obs/obs/gridded/metdata/metdata_tmax_CONUS_monthly.mat');
tmax=tmaxm.data(:,:,:,year-1978);
load /data/obs/obs/gridded/prism/lt71m/800m_grids/deltaoffsets_19792008 deltatmax
for j=1:12
tmaxdata(:,:,j)=interp2(lona,lata,tmax(:,:,j),lon,lat);
end
clear tmax
data=tmaxdata+deltatmax-273.15;
end

case 4,
if year<2010
tmin=matfile(['/data/obs/obs/gridded/prism/lt71m/800m_grids/tmin/tmin_800mPRISM_',num2str(year),'_CONUS.mat']);
data=tmin.data;
else
tminm=matfile('/data/obs/obs/gridded/metdata/metdata_tmin_CONUS_monthly.mat');
tmin=tminm.data(:,:,:,year-1978);
load /data/obs/obs/gridded/prism/lt71m/800m_grids/deltaoffsets_19792008 deltatmin
for j=1:12
tmindata(:,:,j)=interp2(lona,lata,tmin(:,:,j),lon,lat);
end
clear tmin
data=tmindata+deltatmin-273.15;
end

case 5,
if year<2010
ppt=matfile(['/data/obs/obs/gridded/prism/lt71m/800m_grids/ppt/ppt_800mPRISM_',num2str(year),'_CONUS.mat']);
data=ppt.data;
else
pptm=matfile('/data/obs/obs/gridded/metdata/metdata_ppt_CONUS_monthly.mat');
ppt=pptm.data(:,:,:,year-1978);
load /data/obs/obs/gridded/prism/lt71m/800m_grids/deltaoffsets_19792008 deltappt
for j=1:12
pptdata(:,:,j)=interp2(lona,lata,ppt(:,:,j),lon,lat);
end
clear ppt
data=pptdata.*deltappt;
end

case 6,
if year<2010
tdmean=matfile(['/data/obs/obs/gridded/prism/lt71m/800m_grids/tdmean/tdmean_800mPRISM_',num2str(year),'_CONUS.mat']);
data=tdmean.data;
else
sphm=matfile('/data/obs/obs/gridded/metdata/metdata_sph_CONUS_monthly.mat');
sph=sphm.data(:,:,:,year-1978);
cd ..
tdmean=dewpointZ(sph,repmat(el2,[1 1 12]));
cd SNOW
clear sph
for j=1:12
tdmeandata(:,:,j)=interp2(lona,lata,tdmean(:,:,j),lon,lat);
end
load /data/obs/obs/gridded/prism/lt71m/800m_grids/deltaoffsets_19792008 deltatdew
data=tdmeandata+deltatdew-273.15;
end

case 7,
sradm=matfile('/data/obs/obs/gridded/metdata/metdata_srad_CONUS_monthly.mat');
srad=sradm.data(:,:,:,year-1978);
sradmax=max(max(sradm.data(:,:,:,1:38),[],2),[],4);
sradmin=min(min(sradm.data(:,:,:,1:38),[],2),[],4);
srad=(srad-repmat(sradmin,[1 1386 1]))./repmat(sradmax-sradmin,[1 1386 1]);
for j=1:12
data(:,:,j)=interp2(lona,lata,srad(:,:,j),lon,lat);
end
data=1-data;
end
data=single(data);
