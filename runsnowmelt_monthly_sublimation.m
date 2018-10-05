function [SnowMelt,SnowWaterEq,SnowDepth,Sublimation,SnowTemp,E,Energy,Albedo]=runsnowmelt(cloudiness,tmax,tmin,precip_mm,solar,tdmean,vs,elev,subornot,startingSnow,startingSnowDepth_m);	


    %inputs: cloudiness fraction, Tmax, Tmin, precip_mm, Solar, tdmean, vs,
	%elev

    Sublimation=0;
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

if nargin==9
    startingSnowDepth_m=0;
    startingSnowDensity_kg_m3=450;
else
 SnowWaterEqYest = startingSnow;
 startingSnowDensity_kg_m3 = startingSnow/startingSnowDepth_m* WaterDens;%
end    


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

    f=find(NewSnow<0);NewSnow(f)=0;
    
for i=1:length(vs) rh(i)=(log((10+0.001)/0.001)*log((2+0.0002)/0.0002))/(0.41*0.41*vs(i)*86400);end
   % cloudiness 		<- EstCloudiness(Tmax_C,Tmin_C)
	
    AE 				= atmosphericemissivity(Tav, cloudiness,tdmean,elev);%	# (-) Atmospheric Emissivity
    
%#  New Variables	:
    SnowTemp 		=zeros(size(precip_mm));;% 		# Degrees C
	rhos 		= calcsatvap2(SnowTemp);%	# 	vapor density at surface (kg/m3)
	rhoa 		= calcsatvap2(tdmean);%		#	vapor density of atmoshpere (kg/m3) 
	DCoef 		=zeros(size(precip_mm));;%				#   Density Coefficient (-) (Simplified version)
	SnowMelt 		=zeros(size(precip_mm));%				#  (m)
	Albedo 			=groundAlbedo*ones(size(precip_mm));;% This will change for days with snow
	TE 				=SurfEmissiv*ones(size(precip_mm));%	(-) Terrestrial Emissivity
	
%##	Energy Terms
	H 		=zeros(size(precip_mm));%	#	Sensible Heat exchanged (kJ/m2/d)
	E 		=zeros(size(precip_mm));%	#	Vapor Energy	(kJ/m2/d)
	S 		=zeros(size(precip_mm));%	#	Solar Radiation (kJ/m2/d)
	Energy 	=zeros(size(precip_mm));%	# Net Energy (kJ/m2/d)
	Lt 		=zeros(size(precip_mm));%	#	Terrestrial Longwave Radiation (kJ/m2/d)
	La 		= longwave(AE, Tav);%					#	Atmospheric Longwave Radiation (kJ/m2/d)
	G 		= 173;%								#	Ground Condution (kJ/m2/d) 
	P 		= Cw .* R_m .* Tav/30;%					# 	Precipitation Heat (kJ/m2/d), divide precipitation by thirty since it already contains amount

    
    
    
%##  Initial Values.  
	SnowWaterEqYest = startingSnowDepth_m * startingSnowDensity_kg_m3 / WaterDens;%		
	SnowDepth(1)	= startingSnowDepth_m;	
	if NewSnow(1)-R_m(1)>0;
        Albedo(1)=0.98-(0.98-0.50)*exp(-40*(NewSnow(1)-R_m(1))/30);
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
    Energy(1)=Energy(1)*30;
	if startingSnowDepth_m+NewSnow>0
        SnowDensity(1)=min(450, (startingSnowDensity_kg_m3*startingSnowDepth_m + NewSnowDensity(1)*NewSnow(1))/(startingSnowDepth_m+NewSnow(1)));
    else SnowDensity(1)=450;end
    
	SnowMelt(1)= max(0,min((SnowWaterEqYest+NewSnowWatEq(1)), (Energy(1)-SnowHeatCap*(SnowWaterEqYest+NewSnowWatEq(1))*WaterDens*(0-SnowTemp(1)))/(LatHeatFreez*WaterDens) ) );
	SnowDepth(1)= max(0,(SnowWaterEqYest + NewSnowWatEq(1)-SnowMelt(1))*WaterDens/SnowDensity(1));
	SnowWaterEq(1)= max(0,SnowWaterEqYest-SnowMelt(1)+NewSnowWatEq(1));	
    if subornot
        [Sublimation(1),SnowWaterEq(1),SnowDepth(1),SnowDensity(1)]=est_sub(E(1),SnowWaterEq(1),SnowDepth(1),SnowDensity(1));
    end
    
    
for i=2:length(vs)
    	if NewSnowWatEq(i)-R_m(i)>0;
         Albedo(i)=0.98-(0.98-Albedo(i-1))*exp(-40*((NewSnow(i)-R_m(i))/30));
        else
            if SnowDepth(i-1)<0.15
                Albedo(i) = max(groundAlbedo,0.5+(groundAlbedo-0.85)/10);
            else
                lastalbedo=Albedo(i-1);
                for jj=1:30
                lastalbedo = 0.35-(0.35-0.98)*exp(-1*(0.177+(log(1e-5+(-0.3+0.98)/(-.3+lastalbedo)))^2.16)^0.46);
                end          
                Albedo(i)=lastalbedo;
            end
        end
        S(i)=solar(i)*(1-Albedo(i-1));
        if SnowDepth(i-1)>0 TE(i)=0.97;end
		if(SnowWaterEq(i-1) > 0 | NewSnowWatEq(i) > 0) 
			DCoef(i) = 6.2;
			if SnowMelt(i-1)-NewSnowWatEq(i-1) <=0 
				lasttemp=SnowTemp(i-1);
                lastswe=SnowDensity(i-1)*SnowDepth(i-1);
                for jj=1:30
                lasttemp = max(min(0,tmin(i)),min(0,(lasttemp+min(-lasttemp,Energy(i-1)/30/((lastswe+NewSnow(i)*NewSnowDensity(i)/30)*SnowHeatCap*1000)))));
                lastswe=lastswe+NewSnow(i)*NewSnowDensity(i)/30;  % this will overpredict snow since no melt is assumed
                end
                SnowTemp(i)=lasttemp;
                % force snow temp to be < or = to mean
                SnowTemp(i)=min(SnowTemp(i),Tav(i));
            else
                SnowTemp(i)=0;
            end
        end
        Ri(i)=StabilityCorrection(2,0,SnowTemp(i),tmax(i),vs(i));
        rh(i)=rh(i)/Ri(i);
		rhos(i) = calcsatvap2(SnowTemp(i));
        %rhos(i)=rhos(i)*.98;  % adjust for sat vap over ice
		H(i) = 1.29*(Tav(i)-SnowTemp(i))/rh(i);
		E(i) = lambdaV*(rhoa(i)-rhos(i))/rh(i);
		Lt(i) = longwave(TE(i),SnowTemp(i));
		Energy(i) = S(i) + La(i) - Lt(i) + H(i) + E(i) + G + P(i);
        Energy(i)=Energy(i)*30;
		if Energy(i)>0 k = 2; else k =1;end
		if SnowDepth(i-1)+NewSnow(i)>0  
		SnowDensity(i) = min(450,(((SnowDensity(i-1)+k*30*(450-SnowDensity(i-1))*exp(-DCoef(i)))*SnowDepth(i-1)) + NewSnowDensity(i)*NewSnow(i))/(SnowDepth(i-1)+NewSnow(i)));
        else SnowDensity(i)= 450;end
		SnowMelt(i) = max(0,	min( (SnowWaterEq(i-1)+NewSnowWatEq(i)),(Energy(i)-SnowHeatCap*(SnowWaterEq(i-1)+NewSnowWatEq(i))*WaterDens*(0-SnowTemp(i)))/(LatHeatFreez*WaterDens) )  );
        SnowDepth(i) = max(0,(SnowWaterEq(i-1)+NewSnowWatEq(i)-SnowMelt(i))*WaterDens/SnowDensity(i));%
		SnowWaterEq(i) = max(0,SnowWaterEq(i-1)-SnowMelt(i)+NewSnowWatEq(i));%
        if subornot
            [Sublimation(i),SnowWaterEq(i),SnowDepth(i),SnowDensity(i)]=est_sub(E(i),SnowWaterEq(i),SnowDepth(i),SnowDensity(i));
        end
end
SnowWaterEq=single(SnowWaterEq);SnowMelt=single(SnowMelt);
SnowDepth=single(SnowDepth);    

function [Sublimation,SnowWaterEq,SnowDepth,SnowDensity]=est_sub(EV, SnowWaterEq,SnowDepth,SnowDensity);
    WaterDens =1000;%			# kg/m3
    lambdaS = 2850;%             # (kJ/kg) latent heat of sublimation (this used if SnowTemp<0)

    % sublimation added at end 
    Sublimation=-EV/lambdaS/WaterDens;
    Sublimation=Sublimation*30;  % 30 days in a month
    % modifies SnowDepth and SWE
    
    if SnowWaterEq>Sublimation 
        if SnowWaterEq>0 % no limit
        SnowWaterEq=SnowWaterEq-Sublimation;
        SnowDepth=SnowWaterEq/SnowDensity*WaterDens;
        else
            Sublimation=0;SnowWaterEq=0;SnowDepth=0;SnowDensity=450;
        end
    else
        if SnowWaterEq>0
            Sublimation=SnowWaterEq;
            SnowWaterEq=0;SnowDepth=0;SnowDensity=450;
        else
            Sublimation=0;SnowWaterEq=0;SnowDepth=0;SnowDensity=450;
        end
    end
    