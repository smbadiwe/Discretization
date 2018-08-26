function [cpp1]=cp1(A,B,C,D,E,T,nz,ns)
ns=7
nz=15
Tw=773;
T=Tw;

A(1:ns)=[9.656      29.526    29.342      33.933    25.399      34.942     27.437];
B(1:ns)=[4.283E-01 -8.99E-03 -3.54E-03   -8.42E-03  2.0178E-02 -3.996E-02  4.23E-02];
C(1:ns)=[-2.68E-04  3.81E-05  1.0076E-05  2.99E-05 -3.85E-05    1.92E-04  -1.96E-05];
D(1:ns)=[3.18E-08  -3.26E-08 -4.31E-09   -1.78E-08  3.19E-08   -1.53E-07   3.997E-09];
E(1:ns)=[2.77E-11   8.86E-12  2.59E-13    3.69E-12 -8.76E-12    3.93E-11  -2.99E-13];

%Units for Cpg1 J/mol.K
cpp1(1:nz,1:ns)=0;
for is=1:ns
    cpp1=A+(B.*T)+(C*(T.^2))+(D*(T.^3))+(E*(T.^4))
end
 end