function RR = Rxnrate(P,Tw,yCG_in,yO2_in,yH2O_in,nz)
% yN2,yCH4,yH2,yCO2
R = 8.314; nz=15; P(1:nz)=101325;Tw(1:nz)=773;
% if yCG==0.1
%     yCG = 0.1;
% end
% if yO2 == 0.06;
% %     yO2 = 0.06;
% end
% if yN2 == 0.23;
%     yN2 = 0.23;
% end
% if yH2 == 0
%     yH2 = 1E-03;
% end
% if yCH4 == 0
%     yCH4 = 1E-03;
% end
% if yCO2== 0
%     yCO2 = 1E-03;
% end
% if yH2O==0.61
% yH2O = 1-(ycG+yO2+yN2+yCH4+yH2+yCO2);
% end
yCG_in = 0.1;  yO2_in = 0.06; yN2_in = 0.23; yH2_in = 0.0; yCH4_in = 0.0; yCO2_in = 0.0;
yH2O_in= 1-(yCG_in+yO2_in+yN2_in+yH2_in+yCH4_in+yCO2_in);
%Partial pressures
PCG=P.*yCG_in ;PO2=P.*yO2_in; PN2=P.*yN2_in;  PCO2=P.*yCO2_in; PH2=P.*yH2_in ;PCH4=P.*yCH4_in; PH2O=P-(PCG+PO2+PN2+PCO2+PH2+PCH4)
%Reaction order wrt crude glycerol[-]
m=1.06;
%Reaction order wrt water[-]
n=0.56;
%Reaction order wrt oxygen[-]
o=1.86;
%frequency factor
Apr=(((2.135*10^11))/(60*((101325)^(m+n+o))))% units are [kmolC/kgcat.s.Pa^(m+n+o)
%Activation Energy
An=-93540;%in J/mol
%Rate of reaction[molC/gcat.s] same as [kmolC/kgcat.s]
RR(1:nz)=0;
for i=1:nz
RR(i)=(((Apr)*(exp(An./(R*(Tw(i)))))).*((PCG(i).^m).*((PH2O(i)).^n).*((PO2(i)).^o)))
end
end