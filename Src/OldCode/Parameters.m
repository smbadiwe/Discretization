%Prop.
L=0.045;% in m
% diameter of pellets in m
dp=0.0008; 
%diameter of reactor in m
d=0.0127;
ns=7;%number of species
v(1:ns)=[1 0.65 0 0.7 3.2 0.5 2]; % stoichiometric coefficient[-]
nz=15; %number of cells in reactor
nt=1000; %numebr of time steps
dz=L/nz;
t=15; %units in s for the total time
dt=t/nt;
FA0= 2.10049E-06 ;%in m3/s from the excel graph Mol Avg
Area= pi*(d/4)^2;
ep = 0.38 + 0.073 *(1+(((d/dp-2)^2)/((d/dp)^2)));
T0=298;%units in K
uz = FA0/Area;
Rop=53.5;%8.77;%density of pellet in kg/m3
cpp=379.64; % of catalyst pellet[J/KgK]
pW=1.1e5;%[Pa}
pE=1.01325e5;%[Pa}
Press(1:nz)=pE
Tw=773;%units in K
T(1:nz)=Tw;%units in K
uv=0.204;% in Pa.S
kt=0.062;%W/mKs


z=L
z(1)=0;
for i=2:nz
    z(2)=dz/2
    z(3:nz+1)=(z(2)+(1:nz-1)*dz)
    z(nz+2)=z(nz+1)+dz/2
    zEdge(1:nz+1)=(0:nz)*dz
end
