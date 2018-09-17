function [SnowMelt,SnowWaterEq]=runsnowmelt(cloudiness,tmax,tmin,precip_mm,solar,tdmean,vs,elev);	
    %inputs: cloudiness fraction, Tmax, Tmin, precip_mm, Solar, tdmean, vs,
	%elev


    WaterDens =1000;%			# kg/m3
	lambda = 3.35E5;%			# latent heat of fusion (kJ/m3)
	lambdaV = 2500;%				# (kJ/kg) latent heat of vaporization
	SnowHeatCap = 2.1;%			# kJ/kg/C
	LatHeatFreez = 333.3;%		# kJ/kg
	Cw = 4.2E3;%				# Heat Capacity of Water (kJ/m3/C)
    windHt=10;
    tempHt=2;
	groundAlbedo=0.25;
    SurfEmissiv=0.95;
    startingSnowDepth_m=0;
    startingSnowDensity_kg_m3=450;
    
%##	Converted Inputs :
	Tav = (tmax+tmin)/2;%		# degrees C
	precip_m = precip_mm*0.001;%	 	# precip in m 
	R_m = precip_m;%				# (m) depth of rain
    passnow=runsnow(Tav+273.15,ones(size(Tav)));
    R_m=R_m.*(1-passnow);
    NewSnowWatEq = precip_m-R_m;
	NewSnowDensity = 50+3.4*(Tav+15);%		# kg/m3
	f=find(NewSnowDensity<50);
	NewSnowDensity(f)=50;
    NewSnow = NewSnowWatEq.*WaterDens./NewSnowDensity;%		# m

for i=1:365 rh(i)=(log((10+0.001)/0.001)*log((2+0.0002)/0.0002))/(0.41*0.41*vs(i)*86400);end
   % cloudiness 		<- EstCloudiness(Tmax_C,Tmin_C)
	
    AE 				= atmosphericemissivity(Tav, cloudiness,tdmean,elev);%	# (-) Atmospheric Emissivity
    
%#  New Variables	:
    SnowTemp 		=zeros(size(precip_mm));;% 		# Degrees C
	rhos 		= calcsatvap(SnowTemp);%	# 	vapor density at surface (kg/m3)
	rhoa 		= calcsatvap(tmin);%		#	vapor density of atmoshpere (kg/m3) 
	DCoef 		=zeros(size(precip_mm));;%				#   Density Coefficient (-) (Simplified version)
	SnowDensity 	=450*ones(size(precip_mm));;%			#  (kg/m3)  Max density is 450
	SnowMelt 		=zeros(size(precip_mm));%				#  (m)
	SnowDepth 		=zeros(size(precip_mm));	%	#  (m)
	Albedo 			=groundAlbedo*ones(size(precip_mm));;% This will change for days with snow
	SnowWaterEq 	=zeros(size(precip_mm));%		#  (m) Equiv depth of water
	TE 				=SurfEmissiv*ones(size(precip_mm));%	(-) Terrestrial Emissivity
	
%##	Energy Terms
	H 		=zeros(size(precip_mm));%	#	Sensible Heat exchanged (kJ/m2/d)
	E 		=zeros(size(precip_mm));%	#	Vapor Energy	(kJ/m2/d)
	S 		=zeros(size(precip_mm));%	#	Solar Radiation (kJ/m2/d)
	Energy 	=zeros(size(precip_mm));%	# Net Energy (kJ/m2/d)
	Lt 		=zeros(size(precip_mm));%	#	Terrestrial Longwave Radiation (kJ/m2/d)
	La 		= longwave(AE, Tav);%					#	Atmospheric Longwave Radiation (kJ/m2/d)
	G 		= 173;%								#	Ground Condution (kJ/m2/d) 
	P 		= Cw .* R_m .* Tav;%					# 	Precipitation Heat (kJ/m2/d)

    
%##  Initial Values.  
	SnowWaterEqYest = startingSnowDepth_m * startingSnowDensity_kg_m3 / WaterDens;%		
	SnowDepth(1)	= startingSnowDepth_m;	
	if NewSnow(1)>0;
        Albedo(1)=0.98-(0.98-0.50)*exp(-4*NewSnow(1)*10);
    else
%    f=find(startingSnowDepth_m == 0 & NewSnow==0);
        Albedo(1)=max(groundAlbedo,0.5+(groundAlbedo-0.85)/10);
    end
	S(1)=solar(1)*(1-Albedo(1));
	H(1)= 1.29*(Tav(1)-SnowTemp(1))./rh(1);% 
	E(1) = lambdaV*(rhoa(1)-rhos(1))./rh(1);
	if startingSnowDepth_m>0
    	TE(1)=0.97;
    else TE(1)=SurfEmissiv;end
	Lt(1) = longwave(TE(1),SnowTemp(1));
	Energy(1) = S(1) + La(1) - Lt(1) + H(1) + E(1) + G + P(1);
    
	if startingSnowDepth_m+NewSnow>0
        SnowDensity(1)=min(450, (startingSnowDensity_kg_m3*startingSnowDepth_m + NewSnowDensity(1)*NewSnow(1))/(startingSnowDepth_m+NewSnow(1)));
    else SnowDensity(1)=450;end
    
	SnowMelt(1)= max(0,min((SnowWaterEqYest+NewSnowWatEq(1)), (Energy(1)-SnowHeatCap*(SnowWaterEqYest+NewSnowWatEq(1))*WaterDens*(0-SnowTemp(1)))/(LatHeatFreez*WaterDens) ) );
	SnowDepth(1)= max(0,(SnowWaterEqYest + NewSnowWatEq(1)-SnowMelt(1))*WaterDens/SnowDensity(1));
	SnowWaterEq(1)= max(0,SnowWaterEqYest-SnowMelt(1)+NewSnowWatEq(1));	
	
    
for i=2:365
    	if NewSnow(i)>0;
         Albedo(i)=0.98-(0.98-Albedo(i-1))*exp(-4*NewSnow(1)*10);
        else
            if SnowDepth(i-1)<0.1
                Albedo(i) = max(groundAlbedo,0.5+(groundAlbedo-0.85)/10);
            else
                Albedo(i) = 0.35-(0.35-0.98)*exp(-1*(0.177+(log((-0.3+0.98)/(Albedo(i-1)-0.3)))^2.16)^0.46);
            end
        end
        S(i)=solar(i)*(1-Albedo(i-1));
        if SnowDepth(i-1)>0 TE(i)=0.97;end
		if(SnowWaterEq(i-1) > 0 | NewSnowWatEq(i) > 0) 
			DCoef(i) = 6.2;
			if SnowMelt(i-1) == 0 
				SnowTemp(i) = max(min(0,tmin(i)),min(0,(SnowTemp(i-1)+min(-SnowTemp(i-1),Energy(i-1)/((SnowDensity(i-1)*SnowDepth(i-1)+NewSnow(i)*NewSnowDensity(i))*SnowHeatCap*1000)))));
            end
        end
		rhos(i) = calcsatvap(SnowTemp(i));
		H(i) = 1.29*(Tav(i)-SnowTemp(i))/rh(i);
		E(i) = lambdaV*(rhoa(i)-rhos(i))/rh(i);
		Lt(i) = longwave(TE(i),SnowTemp(i));
		Energy(i) = S(i) + La(i) - Lt(i) + H(i) + E(i) + G + P(i);

		if Energy(i)>0 k = 2; else k =1;end
		if SnowDepth(i-1)+NewSnow(i)>0
		SnowDensity(i) = min(450,((SnowDensity(i-1)+k*30*(450-SnowDensity(i-1))*exp(-DCoef(i)))*SnowDepth(i-1) + NewSnowDensity(i)*NewSnow(i))/(SnowDepth(i-1)+NewSnow(i)));
        else SnowDensity(i)= 450;end
		SnowMelt(i) = max(0,	min( (SnowWaterEq(i-1)+NewSnowWatEq(i)),(Energy(i)-SnowHeatCap*(SnowWaterEq(i-1)+NewSnowWatEq(i))*WaterDens*(0-SnowTemp(i)))/(LatHeatFreez*WaterDens) )  );
		SnowDepth(i) = max(0,(SnowWaterEq(i-1)+NewSnowWatEq(i)-SnowMelt(i))*WaterDens/SnowDensity(i));%
		SnowWaterEq(i) = max(0,SnowWaterEq(i-1)-SnowMelt(i)+NewSnowWatEq(i));%
end
i

    