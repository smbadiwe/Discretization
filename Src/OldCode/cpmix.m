
function [cpg] = cpmix(ns,yCG_in,yO2_in,yH2O_in,yN2_in,yCH4_in,yH2_in,yCO2_in,T,nz)
nz=15
ns=7
Tw=773;
T=Tw
A(1:ns)=[9.656      29.526    29.342      33.933    25.399      34.942     27.437];
B(1:ns)=[4.283E-01 -8.99E-03 -3.54E-03   -8.42E-03  2.0178E-02 -3.996E-02  4.23E-02];
C(1:ns)=[-2.68E-04  3.81E-05  1.0076E-05  2.99E-05 -3.85E-05    1.92E-04  -1.96E-05];
D(1:ns)=[3.18E-08  -3.26E-08 -4.31E-09   -1.78E-08  3.19E-08   -1.53E-07   3.997E-09];
E(1:ns)=[2.77E-11   8.86E-12  2.59E-13    3.69E-12 -8.76E-12    3.93E-11  -2.99E-13];
yCG_in = 0.1;  yO2_in = 0.06; yN2_in = 0.23; yH2_in = 0.0; yCH4_in = 0.0; yCO2_in = 0.0;yH2O_in= 1-(yCG_in+yO2_in+yN2_in+yH2_in+yCH4_in+yCO2_in);
%   yCG(1:nz)=yCG_in; yO2(1:nz)=yO2_in; yN2(1:nz)=yN2_in; yH2O(1:nz)=yH2O_in; yH2(1:nz)=yH2_in; yCH4(1:nz)=yCH4_in; yCO2(1:nz)=yCO2_in;

  %Average molecular weight in[g/mol] or [kg/kmol]
MCG=69; MO2=32; MN2=28; MH2O=18;  MH2=2;  MCH4=16;  MCO2=44;%Molecular weight of species Kg/kmol
% for i=1:nz
%  Mt=MCG.*yCG(i) +MO2.*yO2(i) +MN2.*yN2(i) +MH2O.*yH2O(i) +MCH4.*yCH4(i) +MH2.*yH2(i) +MCO2.*yCO2(i)
% end
Mt=MCG*yCG_in +MO2*yO2_in +MN2*yN2_in +MH2O*yH2O_in +MCH4*yCH4_in +MH2*yH2_in +MCO2*yCO2_in
  %specific heat capacity of the mixture is given as

 %   cp(1:nz)=(abs(((yCO2*cpg1(7)1)+(yH2*cpg1(6))+(yH2*cpg1(5)))-(yH2O*cpg1(4))-(yN2*cpg1(3))-(yO2*cpg1(2))-(yCG*cpg1(1)))) % units in J/molK
cpg(1:nz,1:ns)=0;
 for is=1:ns
    cpg1=cp1(A,B,C,D,E,T)
end
cp=(abs(((yCO2_in*cpg1(7))+(yH2_in*cpg1(6))+(yH2_in*cpg1(5)))-(yH2O_in*cpg1(4))-(yN2_in*cpg1(3))-(yO2_in*cpg1(2))-(yCG_in*cpg1(1)))) % units in J/molK
  cpg=((cp*1000)./Mt)%units in J/KgK

%To convert cpg(J/molK) to KJ/KgK
%(J/molK)*(mol/g)*(1000g/1kg)*(1KJ/1000J)
end