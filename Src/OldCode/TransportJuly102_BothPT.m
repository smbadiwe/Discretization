%%%Geometry
function TransportJuly102_BothPT ()
L=0.045;%0.25 % in m
nz=15;
dz=L/nz;
z(1)=0;
z(2)=dz/2;
z(3:nz+1)=(z(2)+(1:nz-1)*dz);
z(nz+2)=z(nz+1)+dz/2;
zEdge(1:nz+1)=(0:nz)*dz;
% of catalyst pellet[J/KgK]
cpp=379.64; 

%The number of species
ns=7;

%ns=[cG O2 N2 H2O H2 CH4 CO2]
%Constant for specific heat capacity[J/molK]
A(1:ns)=[9.656      29.526    29.342      33.933    25.399      34.942     27.437];
B(1:ns)=[4.283E-01 -8.99E-03 -3.54E-03   -8.42E-03  2.0178E-02 -3.996E-02  4.23E-02];
C(1:ns)=[-2.68E-04  3.81E-05  1.0076E-05  2.99E-05 -3.85E-05    1.92E-04  -1.96E-05];
D(1:ns)=[3.18E-08  -3.26E-08 -4.31E-09   -1.78E-08  3.19E-08   -1.53E-07   3.997E-09];
E(1:ns)=[2.77E-11   8.86E-12  2.59E-13    3.69E-12 -8.76E-12    3.93E-11  -2.99E-13];

% %Constant for viscosity calculations in uP
% aa(1:ns)=[-23.119  44.224  42.606  -36.826  27.758  3.844  11.811];
% b(1:ns)=[2.89E-01  5.62E-01  4.75E-01  4.29E-01  2.12E-01  4.011E-01  4.984E-01];
% c(1:ns)=[-3.428E-05 -1.13E-04  -9.88E-05  -1.62E-05  -3.28E-05  -1.43E-04  -1.085E-04];

% stoichiometric coefficient[-]
v(1:ns)=[1 0.65 0 0.7 3.2 0.5 2];

%standard heat of formation of species[kJ/mol]
Hf(1:ns)=[-582.8 0 0 -241.8 0 -74.85 -393.51];

%Standard heat of reaction[KJ/mol]
Hrxno=((v(5)*Hf(5))+(v(6)*Hf(6))+(v(7)*Hf(7)))-(v(1)*Hf(1))-(v(2)*Hf(2))-(v(3)*Hf(3))-(v(4)*Hf(4))

% diameter of pellets in m
dp=0.0008; 

%Molecular weight of species Kg/kmol
M(1:ns)=[69 32 28 18 2 16 44];

%Boundary & Initial Conditions For Species
%Mole fraction of species[mol/mol]
yw(1) = 0.1; 
yw(2) = 0.06;
yw(3) = 0.23;
yw(5) = 0.00001;
yw(6) = 0.00001;
yw(7) = 0.00001;
yw(4) = 1-(yw(1)+yw(2)+yw(3)+yw(5)+yw(6)+yw(7))
y(1:nz,1:ns)=0;
for is=1:ns
y(1:nz,is)=yw(is)
end

%Average molecular weight of the mixture in Kg/kmol
Mavg=M(1)*yw(1)+M(2)*yw(2)+M(3)*yw(3)+M(4)*yw(4)+M(5)*yw(5)+M(6)*yw(6)+M(7)*yw(7)

W(1:nz,1:ns)=0;
WW(1:ns)=0;
%Mass fraction of species
for is=1:ns
WW(is)=(M(is)*yw(is))./Mavg'
end
W(1:nz,1:ns)=WW(is);

%diameter of reactor in m
d=0.0127;
%universal gas constant in J/kmolK
R=8.314*1000; 
Ru=8.314; %[J/molK]
%inlet temeperature in K
Tambient=298;
Tw=773; %in K
T(1:nz)=773;
T0=Tambient;
% T0=500;

%density of pellet in kg/m3
Rop=8.77;%8.77;%53.5;
%effective thermal conductivity in W/mK where k(1) is at the edge
k=0.062;
%heat transfer coefficent in W/m2K
h=14;

%Boundary & Initial Conditions For Flow and Heat
pW=1.1e5;%[in Pa]
pE=1.01325e5;%[Pa}
Press(1:nz)=pE;
%Time
t=15;%TotalTime in s
dt=1E-04;
nt=15 %150000;

%Properties
%porosity
ep=0.38+0.073*(1+((((d/dp)-2)^2)/((d/dp)^2)))
%density of gas mixture in kg/m3
Rog(1:nz)=((Press(1:nz)).*Mavg)/(R*(T(1:nz)))

%for viscosity(microP(10E-06P) formula changed view excel 
% uv1(1:nz,1:ns)=0;
% for is=1:ns
%     uv1(1:nz,is)=(Visc(nz,aa(is),b(is),c(is),T)*1E-06)%[micro-Poise)  
% end
% %Conversion to Poise
% uvmix(1:nz)=(((y(:,7).*uv1(:,7))+(y(:,6).*uv1(:,6))+(y(:,5).*uv1(:,5)))+(y(:,4).*uv1(:,4))+(y(:,3).*uv1(:,3))+(y(:,2).*uv1(:,2))+(y(:,1).*uv1(:,1)))
uv(1:nz)=0.204;%Pa.s %conversion to Pa.s

%units are in J/kmolK for the heat of reaction
 %effectiveness factor for the rate
nk=0.98;
%Reaction order wrt crude glycerol[-]
m=1.06;
%Reaction order wrt water[-]
n=0.56;
%Reaction order wrt oxygen[-]
o=1.86;
%frequency factor  units are [kmolC/kgcat.s.Pa^(m+n+o)
Apr=(((2.135*10^11))/(60*((101325)^(m+n+o))))
% % %Activation Energy[J/mol] this is because R is also in Kmol
An=-93540;

%inlet flowrate
FAO=4.167E-9;

%Area
Area=(3.142*(d^2))/4;
uz_in=FAO/Area
% viscous pressure drop Units [kg/m4s]
kD(1:nz)=(150*uv(1:nz).*((1-ep)^2))/((dp^2)*(ep^3))%229188.56 

%kinetic pressure drop [kg/m4]
kV=(1.75*Rog(1:nz)*(1-ep))./dp*(ep^3) 
%velocity
vel(1)=(-1*(Press(1)-pW))./(kD(1)*(0.5*dz));
vel(2:nz)=(-1*(Press(2:nz)-Press(1:nz-1)))./(kD(2:nz)*dz)
vel(nz+1)=(-1*(pE-Press(nz)))./(kD(nz)*(0.5*dz))

%The average velocity[m/s]
velAvg(1:nz)=(vel(1:nz)+vel(2:nz+1))/2

%sum of viscous and kinetic term where  kp(1) is at the edge-inverse
Kp(1:nz)=1./(kD+(kV.*(uz_in))) %BUG:: Which one is being used to compute Kp?
% Kp(1:nz)=1./(kD+(kV.*(velAvg(1:nz))))
SE(1:nz)=0;
P(1:nz,is)=0;
for is=1:ns
P(1:nz,is)=yw(is).*(Press(1:nz)')
end
%Initialize matrix for Press & Temp in density eqn
Ai(1:nz,1:2*nz)=0.0;
bp(1:nz)=0.0;

%Initialize global matrix matrix for Press & Temp
AG(1:2*nz,1:2*nz)=0.0;
bG(1:2*nz)=0.0;

%for specific heat capacity of gas in J/kgK  where cp(1) is at the edge cpg(1:nz+1)=2020;
H1(1:nz,1:ns)=0;%(J/mol)
cP(1:nz,1:ns)=0;%units in [J/molK] and [J/KgK]
for is=1:ns
[cP(1:nz,is),H1(1:nz,is)]=SpHeat(nz,A(is),B(is),C(is),D(is),E(is),T,T0)
end
%cpmolmix units KJ/KmolK
cPmolmix(1:nz)=(abs(((y(:,7).*cP(:,7))+(y(:,6).*cP(:,6))+(y(:,5).*cP(:,5)))-(y(:,4).*cP(:,4))-(y(:,3).*cP(:,3))-(y(:,2).*cP(:,2))-(y(:,1).*cP(:,1))))
cpg(1:nz)=((cPmolmix(1:nz)*1000)./Mavg)%units in J/KgK

%where H1(1:nz) is the same as integral CpdT
H1_mix(1:nz)=(v(7)*(H1(:,7)))+(v(6)*(H1(:,6)))+((v(5)*H1(:,5)))-(v(4)*(H1(:,4)))-(v(3)*(H1(:,3)))-(v(2)*(H1(:,2)))-(v(1)*(H1(:,1)))
%Heat of reaction term
Hrxn(1:nz)=(Hrxno)+((H1_mix(1:nz))/1000)% in KJ/mol

%Rate of reaction[molC/gcat.s] same as [kmolC/kgcat.s]
RR(1)=(((Apr)*(exp(An./(Ru*(Tw'))))).*((P(1,1)).^m).*((P(1,4)).^n).*((P(1,2)).^o))
RR(2:nz)=0;
RR0=RR(1,1);

%Diffusion coefficient
Dm(1:ns)=[3.8E-5  1.03E-04  0  1.14E-04   3.64E-04  1.15E-04  7.94E-05];

 %Axial dispersion coefficient in m2/s
 Dli(1:nz,is)=0;
 for is=1:ns
 Dli(1:nz,is)=0.73.*(Dm(is))+((0.5.*velAvg(1:nz)*dp)./(1+((9.49.*(Dm(is)))./((velAvg(1:nz)).*dp))))
 end
 
 %Source term for species equation
SW(1:nz,1:ns)=0;
 for is=1:ns
 SW(1:nz,is)=(1-ep)*Rop.*(M(is)).*(nk.*(RR(1:nz)).*(v(is)))
 end
 SW;
%1st term in temperature equation
a=(1-ep)*cpp*Rop;

Ap(1:nz,1:nz)=0.0;
ApT(1:nz,1:nz)=0.0;
AT(1:nz,1:nz)=0.0;
bT(1:nz)=0.0;
bp(1:nz)=0.0;
Aw1(1:nz,1:nz)=0.0;
bW1(1:nz)=0.0;
Aw2(1:nz,1:nz)=0.0;
bW2(1:nz)=0.0;
Aw4(1:nz,1:nz)=0.0;
bW4(1:nz)=0.0;
Aw5(1:nz,1:nz)=0.0;
bW5(1:nz)=0.0;
Aw6(1:nz,1:nz)=0.0;
bW6(1:nz)=0.0;
Aw7(1:nz,1:nz)=0.0;
bW7(1:nz)=0.0;
y(1:nz,is)=0;
%Time Loop

for it=2:nt
    it
    if(it==1)
SE(1:nz)=0;
    elseif(it==it+1)
        SE=(((1-ep)*(Rop)*((-1*Hrxn)).*nk.*(RR(1:nz))))*1E6
    else (it==nt)
        (((1-ep)*(Rop)*((-1*Hrxn)).*nk.*(RR(1:nz))))*1E6
    end
    %Density Equation
    Ap(1,1)=((ep*(Rog(1)))/((Press(1))*(dt)))+(((Kp(1))*(Rog(1)))/(0.5*((dz)^2)))+(((Kp(2))*(Rog(2)))/((dz)^2));
    Ap(1,2)=-(((Kp(2))*(Rog(2)))/((dz)^2));
    bp(1)=(((Kp(1))*(Rog(1)))/(0.5*(dz^2)))*(pW);
    Ap(nz,nz-1)=-(((Kp(nz-1))*(Rog(nz-1)))/((dz)^2));
    Ap(nz,nz)=((ep*(Rog(nz)))/((Press(nz))*(dt)))+(((Kp(nz-1))*(Rog(nz-1)))/((dz)^2))+(((Kp(nz))*(Rog(nz)))/(0.5*(dz)^2));
    bp(nz)=(((Kp(nz))*(Rog(nz))*pE)/(0.5*((dz)^2)));
    for ip=2:nz-1
        Ap(ip,ip)=((ep*(Rog(ip)))/((Press(ip))*dt))+(((Kp(ip-1))*(Rog(ip-1)))/(dz^2))+(((Kp(ip+1))*(Rog(ip+1)))/((dz)^2));
        Ap(ip,ip+1)=-(((Kp(ip+1))*(Rog(ip+1)))/((dz)^2));
        Ap(ip,ip-1)=-(((Kp(ip-1))*(Rog(ip-1)))/((dz)^2));
        bp(ip)=0;
    end
    
    Ct(1:nz)=-(ep*(Rog(1:nz)))./(T(1:nz)*dt);
%     Ct'
  eye(nz);
    Ap1(1:nz,1:nz)=eye(nz);
    for i=1:nz
        ApT(i,i)=Ct(i).*Ap1(i,i);
    end
    Ap
    ApT
    bp
    %=zeros(nz)
% ApT = Ct'.*eye(nz)
    ApTT=[Ap,ApT]
    
   %Energy Equation
    AT(1,1)=((a+(ep*(cpg(1))*(Rog(1))))/(dt))+(((Rog(1))*(cpg(1))*(vel(2)))/(dz))+((3*k)/((dz)^2));
    AT(1,2)=-(k/((dz)^2));
    bT(1)=(SE(1))+(((a+(ep*(cpg(1))*(Rog(1))))/dt)*(T(1)))+((((Rog(1))*(cpg(1))*(vel(1)))/(dz))*(Tw))+((k/(0.5*((dz)^2)))*(Tw));
    AT(nz,nz)=((a+(ep*(cpg(nz))*(Rog(nz))))/(dt))+(((Rog(nz))*(cpg(nz))*(vel(nz+1)))/(dz))+(k/((dz)^2));
    AT(nz,nz-1)=-((((Rog(nz-1))*(cpg(nz-1))*(vel(nz)))/(dz))+(k/((dz)^2)));
    bT(nz)=(SE(nz))+(((a+(ep*(Rog(nz))*(cpg(nz))))/(dt))*(T(nz)));
   
    for iT=2:nz-1
        AT(iT,iT)=((a+(ep*(cpg(iT))*(Rog(iT))))/(dt))+(((Rog(iT))*(cpg(iT))*(vel(iT+1)))/(dz))+((2*k)/((dz)^2));
        AT(iT,iT-1)=-((((Rog(iT-1))*(cpg(iT-1))*(vel(iT)))/(dz))+(k/((dz)^2)));
        AT(iT,iT+1)=-(k/((dz)^2));
        bT(iT)=(SE(iT))+(((a+(ep*(cpg(iT))*(Rog(iT))))/(dt))*(T(iT))) ;
    end
    AT
    Atl=zeros(nz);
    ATT=[Atl,AT]
    AG=[ApTT;ATT]
    bG=[bp,bT]
    
    %Solve global system
    PT=sparse(AG)\sparse(bG');
    Press(1:nz)=PT(1:nz);
    T(1:nz)=PT(nz+1:2*nz);
    PressPlusBoundary=[pW Press pE]'
    TempPlusBoundary=[Tw T T(nz)]'
   
%     %species equation
% Aw1(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(1))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(2,1)))/((dz)^2))+((ep*(Rog(1))*(Dli(1,1)))/(0.5*(dz)^2));
% Aw1(1,2)=-((ep*(Rog(2))*(Dli(2,1)))/((dz)^2));
% bW1(1)=(SW(1,1))+(((ep*Rog(1))/(dt))*(W(1,1)))+(((Rog(1)*vel(1))/(dz))*(WW(1)))+(((ep*Rog(1)*Dli(1,1))/(0.5*(dz)^2))*(WW(1)));
% Aw1(nz,nz-1)=-((((Rog(nz-1))*(vel(nz)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,1)))/((dz)^2)));
% Aw1(nz,nz)=((ep*(Rog(nz)))/(dt))+(((Rog(nz))*(vel(nz+1)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,1)))/((dz)^2));
% bW1(nz)=SW(nz,1)+(((ep*Rog(nz))/(dt))*(W(nz,1)));
% for i=2:nz-1
% Aw1(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(i+1,1)))/((dz)^2))+((ep*(Rog(i-1))*(Dli(i-1,1)))/((dz)^2));
% Aw1(i,i-1)=-((((Rog(i-1))*(vel(i)))/(dz))+((ep*(Rog(i-1))*(Dli(i-1,1)))/((dz)^2)));
% Aw1(i,i+1)=-((ep*(Rog(i+1))*(Dli(i+1,1)))/((dz)^2));
% bW1(i)=SW(i,1)+(((ep*Rog(i))/(dt))*(W(i,1)));
% end
% Aw1
% bW1
% %species 2 equation
% Aw2(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(1))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(2,2)))/((dz)^2))+((ep*(Rog(1))*(Dli(1,2)))/(0.5*(dz)^2));
% Aw2(1,2)=-((ep*(Rog(2))*(Dli(2,2)))/((dz)^2));
% bW2(1)=(SW(1,2))+(((ep*Rog(1))/(dt))*(W(1,2)))+(((Rog(1)*vel(1))/(dz))*(WW(2)))+(((ep*Rog(1)*Dli(1,2))/(0.5*(dz)^2))*(WW(2)));
% Aw1(nz,nz-1)=-((((Rog(nz-1))*(vel(nz)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,2)))/((dz)^2)));
% Aw2(nz,nz)=((ep*(Rog(nz)))/(dt))+(((Rog(nz))*(vel(nz+1)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,2)))/((dz)^2));
% bW2(nz)=SW(nz,2)+(((ep*Rog(nz))/(dt))*(W(nz,2)));
% for i=2:nz-1
% Aw2(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(i+1,2)))/((dz)^2))+((ep*(Rog(i-1))*(Dli(i-1,2)))/((dz)^2));
% Aw2(i,i-1)=-((((Rog(i-1))*(vel(i)))/(dz))+((ep*(Rog(i-1))*(Dli(i-1,2)))/((dz)^2)));
% Aw2(i,i+1)=-((ep*(Rog(i+1))*(Dli(i+1,2)))/((dz)^2));
% bW2(i)=SW(i,2)+(((ep*Rog(i))/(dt))*(W(i,2)));
% end
% Aw2;
% bW2;
% 
%  %for N2
%  W(1:nz,3)=WW(3);
%  
% %species 4 equation 
% Aw4(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(1))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(2,4)))/((dz)^2))+((ep*(Rog(1))*(Dli(1,4)))/(0.5*(dz)^2));
% Aw4(1,2)=-((ep*(Rog(2))*(Dli(2,4)))/((dz)^2));
% bW4(1)=(SW(1,4))+(((ep*Rog(1))/(dt))*(W(1,4)))+(((Rog(1)*vel(1))/(dz))*(WW(4)))+(((ep*Rog(1)*Dli(1,4))/(0.5*(dz)^2))*(WW(4)));
% Aw4(nz,nz-1)=-((((Rog(nz-1))*(vel(nz)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,4)))/((dz)^2)));
% Aw4(nz,nz)=((ep*(Rog(nz)))/(dt))+(((Rog(nz))*(vel(nz+1)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,4)))/((dz)^2));
% bW4(nz)=SW(nz,4)+(((ep*Rog(nz))/(dt))*(W(nz,4)));
% for i=2:nz-1
% Aw4(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(i+1,4)))/((dz)^2))+((ep*(Rog(i-1))*(Dli(i-1,4)))/((dz)^2));
% Aw4(i,i-1)=-((((Rog(i-1))*(vel(i)))/(dz))+((ep*(Rog(i-1))*(Dli(i-1,4)))/((dz)^2)));
% Aw4(i,i+1)=-((ep*(Rog(i+1))*(Dli(i+1,4)))/((dz)^2));
% bW4(i)=SW(i,4)+(((ep*Rog(i))/(dt))*(W(i,4)));
% end
% Aw4;
% bW4;
% 
%  %species 5 equation 
% Aw5(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(1))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(2,5)))/((dz)^2))+((ep*(Rog(1))*(Dli(1,5)))/(0.5*(dz)^2));
% Aw5(1,2)=-((ep*(Rog(2))*(Dli(2,5)))/((dz)^2));
% bW5(1)=(SW(1,5))+(((ep*Rog(1))/(dt))*(W(1,5)))+(((Rog(1)*vel(1))/(dz))*(WW(5)))+(((ep*Rog(1)*Dli(1,5))/(0.5*(dz)^2))*(WW(5)));
% Aw5(nz,nz-1)=-((((Rog(nz-1))*(vel(nz)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,5)))/((dz)^2)));
% Aw5(nz,nz)=((ep*(Rog(nz)))/(dt))+(((Rog(nz))*(vel(nz+1)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,5)))/((dz)^2));
% bW5(nz)=SW(nz,5)+(((ep*Rog(nz))/(dt))*(W(nz,5)));
% for i=2:nz-1
% Aw5(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(i+1,5)))/((dz)^2))+((ep*(Rog(i-1))*(Dli(i-1,5)))/((dz)^2));
% Aw5(i,i-1)=-((((Rog(i-1))*(vel(i)))/(dz))+((ep*(Rog(i-1))*(Dli(i-1,5)))/((dz)^2)));
% Aw5(i,i+1)=-((ep*(Rog(i+1))*(Dli(i+1,5)))/((dz)^2));
% bW5(i)=SW(i,5)+(((ep*Rog(i))/(dt))*(W(i,5)));
% end
% Aw5;
% bW5;
% 
% %species 6 equation 
%  Aw6(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(1))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(2,6)))/((dz)^2))+((ep*(Rog(1))*(Dli(1,6)))/(0.5*(dz)^2));
% Aw6(1,2)=-((ep*(Rog(2))*(Dli(2,6)))/((dz)^2));
% bW6(1)=(SW(1,6))+(((ep*Rog(1))/(dt))*(W(1,6)))+(((Rog(1)*vel(1))/(dz))*(WW(6)))+(((ep*Rog(1)*Dli(1,6))/(0.5*(dz)^2))*(WW(6)));
% Aw6(nz,nz-1)=-((((Rog(nz-1))*(vel(nz)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,6)))/((dz)^2)));
% Aw6(nz,nz)=((ep*(Rog(nz)))/(dt))+(((Rog(nz))*(vel(nz+1)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,6)))/((dz)^2));
% bW6(nz)=SW(nz,6)+(((ep*Rog(nz))/(dt))*(W(nz,6)));
% for i=2:nz-1
% Aw6(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(i+1,6)))/((dz)^2))+((ep*(Rog(i-1))*(Dli(i-1,6)))/((dz)^2));
% Aw6(i,i-1)=-((((Rog(i-1))*(vel(i)))/(dz))+((ep*(Rog(i-1))*(Dli(i-1,6)))/((dz)^2)));
% Aw6(i,i+1)=-((ep*(Rog(i+1))*(Dli(i+1,6)))/((dz)^2));
% bW6(i)=SW(i,6)+(((ep*Rog(i))/(dt))*(W(i,6)));
% end
% Aw6;
% bW6;
% 
% %species 7 equation 
%  Aw7(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(1))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(2,7)))/((dz)^2))+((ep*(Rog(1))*(Dli(1,7)))/(0.5*(dz)^2));
% Aw7(1,2)=-((ep*(Rog(2))*(Dli(2,7)))/((dz)^2));
% bW7(1)=(SW(1,7))+(((ep*Rog(1))/(dt))*(W(1,7)))+(((Rog(1)*vel(1))/(dz))*(WW(7)))+(((ep*Rog(1)*Dli(1,7))/(0.5*(dz)^2))*(WW(7)));
% Aw7(nz,nz-1)=-((((Rog(nz-1))*(vel(nz)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,7)))/((dz)^2)));
% Aw7(nz,nz)=((ep*(Rog(nz)))/(dt))+(((Rog(nz))*(vel(nz+1)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,7)))/((dz)^2));
% bW7(nz)=SW(nz,7)+(((ep*Rog(nz))/(dt))*(W(nz,7)));
% for i=2:nz-1
% Aw7(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(i+1,7)))/((dz)^2))+((ep*(Rog(i-1))*(Dli(i-1,7)))/((dz)^2));
% Aw7(i,i-1)=-((((Rog(i-1))*(vel(i)))/(dz))+((ep*(Rog(i-1))*(Dli(i-1,7)))/((dz)^2)));
% Aw7(i,i+1)=-((ep*(Rog(i+1))*(Dli(i+1,7)))/((dz)^2));
% bW7(i)=SW(i,7)+(((ep*Rog(i))/(dt))*(W(i,7)));
% end
% Aw7;
% bW7;
% 
%    %Solve the mass fraction global system
%     W(1:nz,1)=sparse(Aw1)\sparse(bW1')
%     Species1PlusBoundary=[WW(1) W(1:nz,1)' W(nz,1)]';
%      W(1:nz,2)=sparse(Aw2)\sparse(bW2')
%     Species2PlusBoundary=[WW(2) W(1:nz,2)' W(nz,2)]';
%     Species3PlusBoundary=[WW(3) W(1:nz,3)' W(nz,3)]';
%      W(1:nz,4)=sparse(Aw4)\sparse(bW4')
%     Species4PlusBoundary=[WW(4) W(1:nz,4)' W(nz,4)]';
%      W(1:nz,5)=sparse(Aw5)\sparse(bW5')
%     Species5PlusBoundary=[WW(5) W(1:nz,5)' W(nz,5)]';
%      W(1:nz,6)=sparse(Aw6)\sparse(bW6')
%     Species6PlusBoundary=[WW(6) W(1:nz,6)' W(nz,6)]';
%      W(1:nz,7)=sparse(Aw7)\sparse(bW7')
%     Species7PlusBoundary=[WW(7) W(1:nz,7)' W(nz,7)]';
%   WWT=Species1PlusBoundary+Species2PlusBoundary+Species3PlusBoundary+Species4PlusBoundary+Species5PlusBoundary+Species6PlusBoundary+Species7PlusBoundary
% %mole fraction update
% for is=1:ns
% y(1:nz,is)=sparse(W(1:nz,is)'./(M(is)))/sparse((W(1:nz,1)'./(M(1)))+(W(1:nz,2)'./(M(2)))+(W(1:nz,3)'./(M(3)))+(W(1:nz,4)'./(M(4)))+(W(1:nz,5)'./(M(5)))+(W(1:nz,6)'./(M(6)))+(W(1:nz,7)'./(M(7))))
% end
% y(1:nz,3)=yw(3);
% CG=[yw(1) y(1:nz,1)' y(nz,1)]';
% O2=[yw(2) y(1:nz,2)' y(nz,2)]';
% N2=[yw(3) y(1:nz,3)' y(nz,3)]';
% H2O=[yw(4) y(1:nz,4)' y(nz,4)]';
% CH4=[yw(5) y(1:nz,5)' y(nz,5)]';
% H2=[yw(6) y(1:nz,6)' y(nz,6)]';
% CO2=[yw(7) y(1:nz,7)' y(nz,7)]';
% %Average molecular weight of the mixture in Kg/kmol
% Mavg=M(1)*y(:,1)+M(2)*y(:,2)+M(3)*y(:,3)+M(4)*y(:,4)+M(5)*y(:,5)+M(6)*y(:,6)+M(7)*y(:,7);
% 
%Density Update
Rog=((Press(1:nz)).*Mavg)./(R*(T(1:nz)));

%update of partial pressure
for is=1:ns
P(1:nz,is)=y(1:nz,is).*(Press(1:nz)')
end
% % %viscosity update
% % for is=1:ns
% %     uv1(1:nz,is)=Visc(nz,aa(is),b(is),c(is),T);
% % end
% % %Conversion to Poise
% % uvmix(1:nz)=(((yw(7).*uv1(:,7))+(yw(6).*uv1(:,6))+(yw(5).*uv1(:,5)))-(yw(4).*uv1(:,4))-(yw(3).*uv1(:,3))-(yw(2).*uv1(:,2))-(yw(1).*uv1(:,1)))*1E-03;
% % uv(1:nz)=uvmix*0.1; %conversion to Pa.s
% 
% viscous pressure drop Units [kg/m4s]
kD(1:nz)=(150*uv(1:nz).*((1-ep)^2))/((dp^2)*(ep^3));%229188.56 
%kinetic pressure drop [kg/m4]
kV=(1.75*Rog(1:nz)*(1-ep))./dp*(ep^3); 

 %velocity update
vel(1)=(-Kp(1))*(Press(1)-pW)/(0.5*dz);
vel(2:nz)=-(Kp(2:nz)).*(Press(2:nz)-Press(1:nz-1))/dz;
vel(nz+1)=-(Kp(nz))*(pE-Press(nz))/(0.5*dz);

%The average velocity[m/s]
velAvg(1:nz)=(vel(1:nz)+vel(2:nz+1))/2;
% 
%  kp update
Kp(1:nz)=1./(kD(1:nz)+(kV.*(velAvg(1:nz))));
%Cp and Hrxn update
for is=1:ns
[cP(1:nz,is),H1(1:nz,is)]=SpHeat(nz,A(is),B(is),C(is),D(is),E(is),T,T0);
end
% cp update
cPmolmix(1:nz)=(abs(((y(:,7).*cP(:,7))+(y(:,6).*cP(:,6))+(y(:,5).*cP(:,5)))-(y(:,4).*cP(:,4))-(y(:,3).*cP(:,3))-(y(:,2).*cP(:,2))-(y(:,1).*cP(:,1))));
cpg(1:nz)=((cPmolmix*1000)./Mavg');%units in J/KgK
% %Heat of reaction update
H1_mix(1:nz)=(v(7)*(H1(:,7)))+(v(6)*(H1(:,6)))+((v(5)*H1(:,5)))-(v(4)*(H1(:,4)))-(v(3)*(H1(:,3)))-(v(2)*(H1(:,2)))-(v(1)*(H1(:,1)))
Hrxn(1:nz)=(Hrxno)+((H1_mix(1:nz))/1000)% in KJ/mol

% % Reaction rate update[molC/gcat.s]
RR(1:nz)=((Apr)*(exp(An./(Ru*(T(1:nz)'))))).*((P(1:nz,1)).^m).*((P(1:nz,4)).^n).*((P(1:nz,2)).^o)

%Sourceterm for temperature equation
SE(1:nz)=(((1-ep)*(Rop)*((-1*Hrxn)).*nk.*(RR(1:nz))))*1E6 % in [J/m3s]

% %Source term for species equation for species
%  for is=1:ns
%  SW(1:nz,is)=(1-ep)*Rop.*(M(is)).*(nk.*(RR(1:nz)).*(v(is)));
%  end
% 
% % clear all
% end
%  
figure  (1)
plot(z',PressPlusBoundary)
xlabel('reactor length (m)') 
ylabel('Pressure (Pa)')
axis([0 0.045 1E5 1.12E5])
figure (2)
plot(zEdge',vel)
xlabel('reactor length (m)') 
ylabel('velocity (m/s)') 
axis([0 0.045 0.0 inf])
figure (3)
plot(z',TempPlusBoundary)
xlabel('reactor length (m)') 
ylabel(' Temperature (K)') 
axis([0 0.045 772 inf])
% figure (4)
% % plot(z',Species1PlusBoundary,Species2PlusBoundary,Species3PlusBoundary,Species4PlusBoundary,Species5PlusBoundary,Species6PlusBoundary,Species7PlusBoundary)
% % plot(z',Species1PlusBoundary,'',z',Species2PlusBoundary,'',z',Species3PlusBoundary,'',z',Species4PlusBoundary,'',z',Species5PlusBoundary,'',z',Species6PlusBoundary,'',z',Species7PlusBoundary,'')
% % xlabel('reactor length (m)') 
% % ylabel('mass fraction')
% % legend('Species1PlusBoundary','Species2PlusBoundary','Species3PlusBoundary','Species4PlusBoundary','Species5PlusBoundary','Species6PlusBoundary','Species7PlusBoundary')
% % axis([0 0.045 0 1])
% plot(z',CG,'',z',O2,'',z',N2,'',z',H2O,'',z',H2,'',z',CH4,'',z',CO2,'')
% xlabel('reactor length (m)') 
% ylabel('mole fraction')
% legend('CG','O2','N2','H2O','H2','CH4','CO2')
% axis([0 0.045 0 1])
% 
% figure (5)
% RR1=[RR0 RR(1:nz) RR(nz)]'
% size(RR1)
% plot(z',RR1)
% xlabel('reactor length (m)') 
% ylabel('Reaction rate')
% axis([0 0.045 0 inf])
% 
% figure (6)
% Hrxn1=[Hrxno Hrxn(1:nz) Hrxn(nz)]'
% size(Hrxn1)
% plot(z',Hrxn1)
% xlabel('reactor length (m)') 
% ylabel('Heat of reaction')
% axis([0 0.045 -inf inf])
% 
% figure (7)
% cpg11=[cpg(1,1) cpg(1:nz) cpg(nz)]'
% size(cpg11)
% plot(z',cpg11)
% xlabel('reactor length (m)') 
% ylabel('specific heat capacity')
% axis([0 0.045 -inf inf])
end
end


%Hi is the integral of cpdT
 function [cpg1,Hi]=SpHeat(nz,A,B,C,D,E,T,T0)
    Hi(1:nz)=0;
    Tao(1:nz)=0;
    T0;
    %Tao=T./T0
 Tao=T-T0;
    cpg1=A+(B.*T)+(C*(T.^2))+(D*(T.^3))+(E*(T.^4));
    Hi=A*Tao+((B/2)*(Tao.^2))+((C/3)*(Tao.^3))+((D/4)*(Tao.^4))+((E/5)*(Tao.^5));
%     Hi=A*T0*(Tao-1)+((B/2)*(T0^2).*(Tao.^2-1))+((C/3)*(T0^3)*(Tao.^3-1))+((D/4)*(T0^4)*(Tao.^4))+((E/5)*(T0^5)*(Tao.^5-1))
 end
%     function[uvn]=Visc(nz,aa,b,c,T)
%     T(1:nz);
%     uvn=(aa)+(b*T)+(c*(T.^2));
%     end
