%%%Geometry
function TransportJuly1021_BothPT ()
L=0.045;%0.25 % in m
nz=15;   dz=L/nz;
z(1)=0;    z(2)=dz/2;    z(3:nz+1)=(z(2)+(1:nz-1)*dz);    z(nz+2)=z(nz+1)+dz/2;   zEdge(1:nz+1)=(0:nz)*dz;
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
yw(1) = 0.1;  yw(2) = 0.06; yw(3) = 0.23; yw(5) = 0.00001; yw(6) = 0.00001; yw(7) = 0.00001;
yw(4) = 1-(yw(1)+yw(2)+yw(3)+yw(5)+yw(6)+yw(7))
y(1:nz,1:ns)=0;
for is=1:ns
y(1:nz,is)=yw(is);
end

%Average molecular weight of the mixture in Kg/kmol
for is=1:nz
Mavg=AvgMolarWt(y(:,1),y(:,2),y(:,3),y(:,4),y(:,5),y(:,6),y(:,7));
end

W(1:nz,1:ns)=0; WW(1:ns)=0;
%Mass fraction of species
for is=1:ns
WW(is)=(M(is)*yw(is))./Mavg(1)';
end
for is=1:ns
    W(1:nz,1:ns)=WW(is);
end

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
Rop=53.5;%8.77;%53.5;
%effective thermal conductivity in W/mK where k(1) is at the edge
k=0.062;
##%heat transfer coefficent in W/m2K
h=14;

%Boundary & Initial Conditions For Flow and Heat
pW=1.1E5;%[in Pa]
pE=1.01325E5;%[Pa}
Press(1:nz)=pE;
%Time
tt=15;%TotalTime in s
dt=1E-04;
nt=10000;%150000;

%Properties
%porosity
ep=0.38+0.073*(1+((((d/dp)-2)^2)/((d/dp)^2)));
%viscosity in Pa.S
uv(1:nz)=0.204;%Pa.s %conversion to Pa.s

%density of gas mixture in kg/m3
Rog(1:nz)=((Press(1:nz)).*Mavg)/(R*(T(1:nz)))
Rog(nz+1)=((pE).*Mavg(nz))/(R*(T(nz)))
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
FAO=2.10049E-06;
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
% velAvg(1:nz)=(vel(1:nz)+vel(2:nz+1))/2

%sum of viscous and kinetic term where  kp(1) is at the edge-inverse
Kp(1:nz)=1./(kD) % BUG:: Where is the "+k_v u_z" part of the eqn? (Eq 11.2 in docx)
 SE(1:nz)=0;%10E9;

P(1:nz,is)=0;
for is=1:ns
P(1:nz,is)=yw(is).*(Press(1:nz)');
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
[cP(1:nz,is),H1(1:nz,is)]=SpHeat(nz,A(is),B(is),C(is),D(is),E(is),T,T0);
end

%cpmolmix units KJ/KmolK
cPmolmix=(abs(((y(:,7).*cP(:,7))+(y(:,6).*cP(:,6))+(y(:,5).*cP(:,5)))-(y(:,4).*cP(:,4))-(y(:,3).*cP(:,3))-(y(:,2).*cP(:,2))-(y(:,1).*cP(:,1))))
cpg=((cPmolmix*1000)./Mavg)%units in J/KgK

%where H1(1:nz) is the same as integral CpdT
H1_mix=(v(7)*(H1(:,7)))+(v(6)*(H1(:,6)))+((v(5)*H1(:,5)))-(v(4)*(H1(:,4)))-(v(3)*(H1(:,3)))-(v(2)*(H1(:,2)))-(v(1)*(H1(:,1)))
%Heat of reaction term
Hrxn=(Hrxno)+((H1_mix)/1000)% in KJ/mol

%Rate of reaction[molC/gcat.s] same as [kmolC/kgcat.s]
RR(1:nz)=(((Apr)*(exp(An./(Ru*(Tw'))))).*((P(1,1)).^m).*((P(1,4)).^n).*((P(1,2)).^o))
% RR(2:nz)=0;
RR0=RR(1,1);

%Diffusion coefficient
Dm(1:ns)=[3.8E-5  1.03E-04  0  1.14E-04   3.64E-04  1.15E-04  7.94E-05];

 %Axial dispersion coefficient in m2/s
 Dli(1:nz+1,is)=0;
 for is=1:ns
 Dli(1:nz+1,is)=(0.73.*(Dm(is)))+((0.5*(vel(1:nz+1))*dp)./(1+(9.49.*(Dm(is)))./((vel(1:nz+1)).*dp)));
 end
 
 %Source term for species equation
SW(1:nz,1:ns)=0;
 for is=1:ns
 SW(1:nz,is)=(1-ep)*Rop.*(M(is)).*(nk.*(RR(1:nz)).*(v(is)));
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
    %it
    
    %Density Equation
    Ap(1,1)=((ep*(Rog(1)))/((Press(1))*(dt)))-(((Kp(1))*(Rog(1)))/(0.5*((dz)^2)))-(((Kp(2))*(Rog(2)))/((dz)^2));
    Ap(1,2)=(((Kp(2))*(Rog(1)))/((dz)^2));
    bp(1)=-(((Kp(1))*(Rog(1)))/(0.5*(dz^2)))*(pW);
    Ap(nz,nz-1)=(((Kp(nz-1))*(Rog(nz-1)))/((dz)^2));
    Ap(nz,nz)=((ep*(Rog(nz)))/((Press(nz))*(dt)))-(((Kp(nz-1))*(Rog(nz-1)))/((dz)^2))-(((Kp(nz))*(Rog(nz)))/(0.5*(dz)^2));
    bp(nz)=-(((Kp(nz))*(Rog(nz))*pE)/(0.5*((dz)^2)));
    for ip=2:nz-1
        Ap(ip,ip)=((ep*(Rog(ip)))/((Press(ip))*dt))-(((Kp(ip-1))*(Rog(ip-1)))/(dz^2))-(((Kp(ip+1))*(Rog(ip+1)))/((dz)^2));
        Ap(ip,ip+1)=(((Kp(ip+1))*(Rog(ip+1)))/((dz)^2));
        Ap(ip,ip-1)=(((Kp(ip-1))*(Rog(ip-1)))/((dz)^2));
        bp(ip)=0;
    end
    
    Ct(1:nz)=-(ep*(Rog(1:nz)))./(T(1:nz)*dt);
    Ct'
    % eye(nz);
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
    
%     SE(1:nz)=10E9; %BUG:: Where is SE set? In 126, it's set to 0s.
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
    %AT
    Atl=zeros(nz);
    ATT=[Atl,AT];
    AG=[ApTT;ATT];
    bG=[bp,bT];
end
%Solve global system
    PT=sparse(AG)\sparse(bG')
    Press(1:nz)=PT(1:nz);
    T(1:nz)=PT(nz+1:2*nz);
    PressPlusBoundary=[pW Press pE]'
    TempPlusBoundary=[Tw T T(nz)]'
    
    %update of parameters
 %velocity update
vel(1)=(-Kp(1))*(Press(1)-pW)/(0.5*dz);
vel(2:nz)=-(Kp(2:nz)).*(Press(2:nz)-Press(1:nz-1))/dz;
vel(nz+1)=-(Kp(nz))*(pE-Press(nz))/(0.5*dz);

%Average molecular weight of the mixture in Kg/kmol
for is=1:nz
Mavg=AvgMolarWt(y(:,1),y(:,2),y(:,3),y(:,4),y(:,5),y(:,6),y(:,7));
end

%Density Update
Rog=((Press(1:nz)).*Mavg)./(R*(T(1:nz)));

%Cp and Hrxn update
for is=1:ns
[cP(1:nz,is),H1(1:nz,is)]=SpHeat(nz,A(is),B(is),C(is),D(is),E(is),T,T0);
end
% cp update
cPmolmix=(abs(((y(:,7).*cP(:,7))+(y(:,6).*cP(:,6))+(y(:,5).*cP(:,5)))-(y(:,4).*cP(:,4))-(y(:,3).*cP(:,3))-(y(:,2).*cP(:,2))-(y(:,1).*cP(:,1))));
cpg=((cPmolmix*1000)./Mavg');%units in J/KgK
% %Heat of reaction update
H1_mix=(v(7)*(H1(:,7)))+(v(6)*(H1(:,6)))+((v(5)*H1(:,5)))-(v(4)*(H1(:,4)))-(v(3)*(H1(:,3)))-(v(2)*(H1(:,2)))-(v(1)*(H1(:,1)))
Hrxn=(Hrxno)+((H1_mix)/1000)% in KJ/mol

%update of partial pressure
for is=1:ns
P(:,is)=y(1:nz,is).*(Press(1:nz)');
end

% % Reaction rate update[molC/gcat.s]
RR=((Apr)*(exp(An./(Ru*(T(1:nz)'))))).*((P(1:nz,1)).^m).*((P(1:nz,4)).^n).*((P(1:nz,2)).^o)

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
end

%Hi is the integral of cpdT
 function [cpg1,Hi]=SpHeat(nz,A,B,C,D,E,T,T0)
    Hi(1:nz)=0;
    Tao(1:nz)=0;
    T0;
 Tao=T-T0;
    cpg1=A+(B.*T)+(C*(T.^2))+(D*(T.^3))+(E*(T.^4));
    Hi=A*Tao+((B/2)*(Tao.^2))+((C/3)*(Tao.^3))+((D/4)*(Tao.^4))+((E/5)*(Tao.^5));
 end
 
 function Mavg=AvgMolarWt(yw1,yw2,yw3,yw4,yw5,yw6,yw7)
 ns=7;
 M(1:ns)=[69 32 28 18 2 16 44];
Mavg=M(1)*yw1 +M(2)*yw2 +M(3)*yw3 +M(4)*yw4 +M(5)*yw5 +M(6)*yw6 +M(7)*yw7;
 end
 
 

