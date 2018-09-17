% need to load solar, tmax, tmin, tdew, vs, solar, ppt, elevation
y=1:3105;
xinc=281;

for ii=25:25

x=1+xinc*(ii-1):xinc*ii;

for i=1:39
tmax=load800_partial(3,1978+i,y,x);
tmin=load800_partial(4,1978+i,y,x);
tmean(:,:,:,i)=tmax/2+tmin/2;clear tmax tmin
end

for i=1:39
[ppt(:,:,:,i),lon,lat,elev]=load800_partial(5,1978+i,y,x);
snowmelt(:,:,:,i)=load800_partial(9,1978+i,y,x);
pet(:,:,:,i)=load800_partial(8,1978+i,y,x);
end

pet=real(pet);
snowmelt=real(snowmelt);
cd ..

passnow=runsnow(tmean+273.15,ones(size(tmean)));
rain=ppt.*(1-passnow);
snowmelt=snowmelt*1000;
rain=rain+snowmelt;clear snowmelt ppt

pet=pet.*(1-passnow);
clear tmean passnow

m=matfile('/home/abatz/SOILS/Polaris_Soils');
xx=m.x;yy=m.y;awc=m.awc_top100;
soil=interp2(xx,yy,awc,lon,lat);
f=find(isnan(soil));soil(f)=100;

rain=shiftdim(reshape(shiftdim(rain,2),12,39,3105*xinc),2);
pet=shiftdim(reshape(shiftdim(pet,2),12,39,3105*xinc),2);

pet(:,12,39)=pet(:,11,39);
rain(:,12,39)=rain(:,11,39);
pet(:,1,39)=pet(:,1,38);

[AET,DEF,RUNOFF,SOILS,DR]=simplehydromodelsnow_grid2_advanced(rain,pet,soil);
AET=shiftdim(reshape(shiftdim(AET,1),12,39,3105,xinc),2);
DEF=shiftdim(reshape(shiftdim(DEF,1),12,39,3105,xinc),2);
cd SNOW
save(['wb_800m_',num2str(ii)],'-v7.3','AET','DEF');
clear AET DEF rain pet soil
ii
end
