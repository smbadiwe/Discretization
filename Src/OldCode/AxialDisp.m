function Dli=AxialDisp(nz,ns,Dm)
% DCG_in,DO2, DN2,DH2O, DH2,DCH4,DCO2
%Diffusion coefficient
DCG_in=3.8E-5; DO2_in=1.03E-04 ; DN2_in= 0; DH2O_in= 1.14E-04 ; DH2_in= 3.64E-04; DCH4_in= 1.15E-04 ;DCO2_in= 7.94E-05;
% Dm(1:ns)=[3.8E-5  1.03E-04  0  1.14E-04   3.64E-04  1.15E-04  7.94E-05];
% Dli(1:nz,1:ns)=0;
% for is=1:ns
%  Dli(1:nz,is)=0.73.*(Dm(is))+((0.5.*velAvg(1:nz)*dp)./(1+((9.49.*(Dm(is)))./((velAvg(1:nz)).*dp))))
%  end
Dli_CG(1:nz)=0;Dli_O2(1:nz)=0;Dli_N2(1:nz)=0;Dli_H2O(1:nz)=0;Dli_H2(1:nz)=0;Dli_CH4(1:nz)=0;Dli_CO2(1:nz)=0;
for i=1:nz
Dli_CG(i)=0.73.*(DCG_in)+((0.5.*velAvg(i)*dp)./(1+((9.49.*(DCG_in))./((velAvg(i)).*dp))));
end
for i=1:nz
Dli_O2(i)=0.73.*(DO2_in)+((0.5.*velAvg(i)*dp)./(1+((9.49.*(DO2_in))./((velAvg(i)).*dp))));
end
for i=1:nz
Dli_N2(i)=0.73.*(DN2_in)+((0.5.*velAvg(i)*dp)./(1+((9.49.*(DN2_in))./((velAvg(i)).*dp))));
end
for i=1:nz
Dli_H2O(i)=0.73.*(DH2O_in)+((0.5.*velAvg(i)*dp)./(1+((9.49.*(DH2O_in))./((velAvg(i)).*dp))));
end
for i=1:nz
Dli_H2(i)=0.73.*(DH2_in)+((0.5.*velAvg(i)*dp)./(1+((9.49.*(DH2_in))./((velAvg(i)).*dp))));
end
for i=1:nz
Dli_CO2(i)=0.73.*(DCO2_in)+((0.5.*velAvg(i)*dp)./(1+((9.49.*(DCO2_in))./((velAvg(i)).*dp))));
end
for i=1:nz
Dli_CH4(i)=0.73.*(DCH4_in)+((0.5.*velAvg(i)*dp)./(1+((9.49.*(DCH4_in))./((velAvg(i)).*dp))));
end
end