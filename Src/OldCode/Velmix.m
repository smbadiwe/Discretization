function uz=Velmix(P,uz_in,Kp,kD,kV,nz)
  FAO_in=2.10049E-06; %units in m3/s
  dp=0.0008; %in m
  d=0.0127 ;%in m
  nz=15
  pE=101325
  P(1:nz)=pE
  L=0.045;%units in m
  dz=L/nz;
  T=773;
  Area=pi*(d/4)^2;%0.0001267; %units in m2
  uz_in=FAO_in/Area
  uv=0.204
  ep = 0.38 + 0.073 *(1+(((d/dp-2)^2)/((d/dp)^2)));
  
  yCG_in = 0.1;  yO2_in = 0.06; yN2_in = 0.23; yH2_in = 0.0; yCH4_in = 0.0; yCO2_in = 0.0;
  yH2O_in= 1-(yCG_in+yO2_in+yN2_in+yH2_in+yCH4_in+yCO2_in);
  %viscous pressure drop
  kD=(150*uv.*((1-ep)^2))/((dp^2)*(ep^3))%229188.56 
  Rog=MixDensity(pE,T,yCG_in,yO2_in,yH2O_in,yN2_in,yCH4_in,yH2_in,yCO2_in,nz)
  %kinetic pressure drop [kg/m4]
  kV=(1.75*Rog*(1-ep))./dp*(ep^3)
  %sum of the kD and kV
  Kp=1./(kD+(kV.*uz_in))
  uz(1)=uz_in;%(-Kp*(P(1)-pW))./(0.5*dz);
  uz(nz+1)=(-Kp(nz)*(pE-P(nz)))./(0.5*dz)
  for i=2:nz
    uz(i)=(-Kp(i).*(P(i)-P(i-1)))./(dz)
  end
  
  %The average velocity[m/s]
  velAvg(1:nz)=(uz(1:nz)+uz(2:nz+1))/2

end