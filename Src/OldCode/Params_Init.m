nz=15; ns=7; nt=1000; t=15;L=0.045;
dz=L/nz; dt=t/nt;
pW=1.1e5; pE=1.01325e5; Tw=773; T=Tw; Press(1:nz)=pE; uz_in= 0.0663; kt=0.062;%W/m2K
DCG_in=3.8E-5; DO2_in=1.03E-04 ; DN2_in= 0; DH2O_in= 1.14E-04 ; DH2_in= 3.64E-04; DCH4_in= 1.15E-04 ;DCO2_in= 7.94E-05;
%Mol fraction
yCG_in = 0.1;  yO2_in = 0.06; yN2_in = 0.23; yH2_in = 0.0; yCH4_in = 0.0; yCO2_in = 0.0;
yH2O_in= 1-(yCG_in+yO2_in+yN2_in+yH2_in+yCH4_in+yCO2_in);
yCG(1:nz)=yCG_in; yO2(1:nz)=yO2_in; yN2(1:nz)=yN2_in; yH2O(1:nz)=yH2O_in; yH2(1:nz)=yH2_in; yCH4(1:nz)=yCH4_in; yCO2(1:nz)=yCO2_in;
uv=0.204;
dp=0.0008;
d=0.0127;
ep = 0.38 + 0.073 *(1+(((d/dp-2)^2)/((d/dp)^2)));
Area=pi*(d/4)^2;
kD=(150*uv.*((1-ep)^2))/((dp^2)*(ep^3))%229188.56 
Rog=MixDensity(Press,T,yCG,yO2,yCO2,yN2,yH2O,yH2,yCH4,nz)
% Rog=MixDensity(pE,T,yCG_in,yO2_in,yH2O_in,yN2_in,yCH4_in,yH2_in,yCO2_in,nz)
%kinetic pressure drop [kg/m4]
kV=(1.75*Rog*(1-ep))./dp*(ep^3)
%sum of the kD and kV
Kp=1./(kD+(kV.*uz_in))

%Mass fraction
WCG_in=0.263; WO2_in = 0.0732; WN2_in = 0.2454; WH2O_in=0.4184; WH2_in=0.0; WCH4_in=0.0; WCO2_in=0.0;
W_in=WCG_in+ WO2_in+ WN2_in+ WH2O_in+ WH2_in+ WCH4_in+ WCO2_in;
WCG(1:nz)=WCG_in; WO2(1:nz)=WO2_in;WN2(1:nz)=WN2_in;WH2O(1:nz)=WH2O_in;WH2(1:nz)=WH2_in;WCH4(1:nz)=WCH4_in;WCO2(1:nz)=WCO2_in;
%Initial transport properties
cpg =cpmix(yCG,yO2,yH2O,yN2,yCH4,yH2,yCO2)
% Rog(1:nz)=0; Dli(1:nz)=0;

vel=Velmix(Press,uz_in,Kp,kD,kV,nz)
velAvg(1:nz)=(vel(1:nz)+vel(2:nz+1))/2
% for is=1:ns
% Dli(1:nz,is)=AxialDisp(nz,ns,DCG,DO2, DN2,DH2O,DH2,DCH4,DCO2)
% end
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

%Average Mass Density of Mixture
MCG=69; MO2=32; MN2=28; MH2O=18;  MH2=2;  MCH4=16;  MCO2=44;%Molecular weight of species Kg/kmol

Mavg=MCG*yCG_in +MO2*yO2_in + MN2*yN2_in +MH2O*yH2O_in +MCH4*yCH4_in +MH2*yH2_in +MCO2*yCO2_in
 yCG(1:nz)=0; yO2(1:nz)=0; yN2(1:nz)=0; yH2O(1:nz)=0; yH2(1:nz)=0; yCH4(1:nz)=0; yCO2(1:nz)=0;
 
 %Mass fraction Conversion to Mole fraction
 for i=1:nz
 yCG(i)=(WCG(i).*MCG)/Mavg
 end
 for i=1:nz
 yO2(i)=(WO2(i).*MO2)/Mavg
 end
 for i=1:nz
 yN2(i)=(WN2(i).*MN2)/Mavg
 end
 for i=1:nz
 yH2O(i)=(WH2O(i).*MH2O)/Mavg
 end
 for i=1:nz
 yH2(i)=(WH2(i).*MH2)/Mavg
 end
 for i=1:nz
 yCH4(i)=(WCH4(i).*MCH4)/Mavg
 end
 for i=1:nz
 yCO2(i)=(WCO2(i).*MCO2)/Mavg
 end
 for i=1:nz
 yCG(i)=(WCG(i).*MCG)/Mavg
 end
y_in=yCG_in+ yO2_in+ yN2_in+ yH2O_in+ yH2_in+ yCH4_in+ yCO2_in;
 
MatrixPress=sparse(ones(nz,nt))*pE
MatrixTemp=sparse(ones(nz,nt))*Tw


SpeciesFractionCG=sparse(ones(nz,nt))*yCG_in
SpeciesFractionO2=sparse(ones(nz,nt))*yO2_in
SpeciesFractionN2=sparse(ones(nz,nt))*yN2_in
SpeciesFractionH2O=sparse(ones(nz,nt))*yH2O_in
SpeciesFractionH2=sparse(ones(nz,nt))*yH2_in
SpeciesFractionCH4=sparse(ones(nz,nt))*yCH4_in
SpeciesFractionCO2=sparse(ones(nz,nt))*yCO2_in
SpeciesFractionAll=sparse(ones(nz,nt))*y_in

MassSpeciesFractionCG=sparse(ones(nz,nt))*WCG_in
MassSpeciesFractionO2=sparse(ones(nz,nt))*WO2_in
MassSpeciesFractionN2=sparse(ones(nz,nt))*WN2_in
MassSpeciesFractionH2O=sparse(ones(nz,nt))*WH2O_in
MassSpeciesFractionH2=sparse(ones(nz,nt))*WH2_in
MassSpeciesFractionCH4=sparse(ones(nz,nt))*WCH4_in
MassSpeciesFractionCO2=sparse(ones(nz,nt))*WCO2_in
MassSpeciesFractionAll=sparse(ones(nz,nt))*W_in


MatrixDens = sparse(ones(nz,nt))*Rog(1);
Matrix_uv = sparse(ones(nz,nt))*uv;
Matrix_kt = sparse(ones(nz,nt))*kt;
Matrixcpg = sparse(ones(nz,nt))*cpg;
MatrixVel = sparse(ones(nz,nt))*uz_in;

Rop = 53.5; % kg/m^3 for catalyst bulk density
catdV = Rop * Area * dz; % Amount of catalyst (kg) in each segment of bed.

RR1 = Rxnrate(Press,T,yCG,yO2,yH2O)
MatrixRR1 = ones(nz,nt)*RR1(1);
