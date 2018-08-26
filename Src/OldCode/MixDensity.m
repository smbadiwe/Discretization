% function [Rog]=MixDensity(pE,T,yCG_in,yO2_in,yH2O_in,yN2_in,yCH4_in,yH2_in,yCO2_in,MCG,MO2,MN2,MH2,MCO2,MH2O,MCH4,Molwt,nz)
function [Rog]=MixDensity(pE,T,yCG_in,yO2_in,yH2O_in,yN2_in,yCH4_in,yH2_in,yCO2_in,nz)
%Average molecular weight in[g/mol] or [kg/kmol]
R=8314
nz=15
pE(1:nz)=101325; Tw=773;
T(1:nz)=Tw

yCG_in = 0.1;  yO2_in = 0.06; yN2_in = 0.23; yH2_in = 0.0; yCH4_in = 0.0; yCO2_in = 0.0;
yH2O_in= 1-(yCG_in+yO2_in+yN2_in+yH2_in+yCH4_in+yCO2_in);

MCG=69; MO2=32; MN2=28; MH2O=18;  MH2=2;  MCH4=16;  MCO2=44;%Molecular weight of species Kg/kmol

Molwt=MCG*yCG_in + MO2*yO2_in + MN2*yN2_in + MH2O*yH2O_in + MCH4*yCH4_in + MH2*yH2_in + MCO2*yCO2_in
% Mt=MCG.*yCG' + MO2.*yO2' + MN2.*yN2'+ MH2O.*yH2O'+ MCH4.*yCH4' + MH2.*yH2' + MCO2.*yCO2'

% Rog=(P*Mt)/(R*T)
% Rog =(pE.*Molwt)./(8314*T);
Rog(1:nz)=0;
 for i=1:nz
Rog(i)=(pE*Molwt)/(R.*T)
 end
     end
