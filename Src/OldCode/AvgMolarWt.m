function Mavg=AvgMolarWt(yCG,yO2,yN2,yH2,yH2O,yCO2,yCH4)
yCG_in = 0.1;  yO2_in = 0.06; yN2_in = 0.23; yH2_in = 0.0; yCH4_in = 0.0; yCO2_in = 0.0;
yH2O_in= 1-(yCG_in+yO2_in+yN2_in+yH2_in+yCH4_in+yCO2_in);

MCG=69; MO2=32; MN2=28; MH2O=18;  MH2=2;  MCH4=16;  MCO2=44;%Molecular weight of species Kg/kmol

Mavg=MCG*yCG_in +MO2*yO2_in +MN2*yN2_in +MH2O*yH2O_in +MCH4*yCH4_in +MH2*yH2_in +MCO2*yCO2_in;
end