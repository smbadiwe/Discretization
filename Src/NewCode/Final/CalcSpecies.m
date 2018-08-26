%species equation 1
Aw1(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(1))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(1,2)))/((dz)^2))+((ep*(Rog(1))*(Dli(1,1)))/(0.5*(dz)^2));
Aw1(1,2)=-((ep*(Rog(2))*(Dli(1,2)))/((dz)^2));
bW1(1)=(SW(1,1))+((ep*(Rog(1))*(W(1,1)))/(dt))+(((Rog(1)*vel(1))/(dz))*(WW(1)))+(((ep*Rog(1)*Dli(1,1))/(0.5*(dz)^2))*(WW(1)));
Aw1(nz,nz-1)=-((((Rog(nz))*(vel(nz)))/(dz))+((ep*(Rog(nz))*(Dli(1,nz)))/((dz)^2)));
Aw1(nz,nz)=((ep*(Rog(nz)))/(dt))+(((RogE)*(vel(nz+1)))/(dz))+((ep*(Rog(nz))*(Dli(1,nz)))/((dz)^2));
bW1(nz)=SW(1,nz)+(((ep*Rog(nz))/(dt))*(W(1,nz)));

for i=2:nz-1
    Aw1(i,i-1)=-((((Rog(i))*(vel(i)))/(dz))+((ep*(Rog(i))*(Dli(1,i)))/((dz)^2)));
    Aw1(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i+1))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(1,i+1)))/((dz)^2))+((ep*(Rog(i))*(Dli(1,i)))/((dz)^2));
    Aw1(i,i+1)=-((ep*(Rog(i+1))*(Dli(1,i+1)))/((dz)^2));
    bW1(i)=SW(1,i)+(((ep*Rog(i))/(dt))*(W(1,i)));
end

%species equation 2
Aw2(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(2))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(2,2)))/((dz)^2))+((ep*(Rog(1))*(Dli(2,1)))/(0.5*(dz)^2));
Aw2(1,2)=-((ep*(Rog(2))*(Dli(2,1)))/((dz)^2));
bW2(1)=(SW(2,1))+((ep*(Rog(1))*(W(2,1)))/(dt))+(((Rog(1)*vel(1))/(dz))*(WW(2)))+(((ep*Rog(1)*Dli(2,1))/(0.5*(dz)^2))*(WW(2)));
Aw2(nz,nz-1)=-((((Rog(nz))*(vel(nz)))/(dz))+((ep*(Rog(nz))*(Dli(2,nz)))/((dz)^2)));
Aw2(nz,nz)=((ep*(Rog(nz)))/(dt))+(((RogE)*(vel(nz+1)))/(dz))+((ep*(Rog(nz))*(Dli(2,nz)))/((dz)^2));
bW2(nz)=SW(2,nz)+(((ep*Rog(nz))/(dt))*(W(2,nz)));

for i=2:nz-1
    Aw2(i,i-1)=-((((Rog(i))*(vel(i)))/(dz))+((ep*(Rog(i))*(Dli(2,i)))/((dz)^2)));
    Aw2(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i+1))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(2,i+1)))/((dz)^2))+((ep*(Rog(i))*(Dli(2,i)))/((dz)^2));
    Aw2(i,i+1)=-((ep*(Rog(i+1))*(Dli(2,i+1)))/((dz)^2));
    bW2(i)=SW(2,i)+(((ep*Rog(i))/(dt))*(W(2,i)));
end

%species equation 3
W(3,1:nz)=WW(3);

%species equation 4
Aw4(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(2))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(4,2)))/((dz)^2))+((ep*(Rog(1))*(Dli(4,1)))/(0.5*(dz)^2));
Aw4(1,2)=-((ep*(Rog(2))*(Dli(4,2)))/((dz)^2));
bW4(1)=(SW(4,1))+((ep*(Rog(1))*(W(4,1)))/(dt))+(((Rog(1)*vel(1))/(dz))*(WW(4)))+(((ep*Rog(1)*Dli(4,1))/(0.5*(dz)^2))*(WW(4)));
Aw4(nz,nz-1)=-((((Rog(nz))*(vel(nz)))/(dz))+((ep*(Rog(nz))*(Dli(4,nz)))/((dz)^2)));
Aw4(nz,nz)=((ep*(Rog(nz)))/(dt))+(((RogE)*(vel(nz+1)))/(dz))+((ep*(Rog(nz))*(Dli(4,nz)))/((dz)^2));
bW4(nz)=SW(4,nz)+(((ep*Rog(nz))/(dt))*(W(4,nz)));

for i=2:nz-1
    Aw4(i,i-1)=-((((Rog(i))*(vel(i)))/(dz))+((ep*(Rog(i))*(Dli(4,i)))/((dz)^2)));
    Aw4(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i+1))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(4,i+1)))/((dz)^2))+((ep*(Rog(i))*(Dli(4,i)))/((dz)^2));
    Aw4(i,i+1)=-((ep*(Rog(i+1))*(Dli(4,i+1)))/((dz)^2));
    bW4(i)=SW(4,i)+(((ep*Rog(i))/(dt))*(W(4,i)));
end

%species equation 5
Aw5(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(2))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(5,2)))/((dz)^2))+((ep*(Rog(1))*(Dli(5,1)))/(0.5*(dz)^2));
Aw5(1,2)=-((ep*(Rog(2))*(Dli(5,2)))/((dz)^2));
bW5(1)=(SW(5,1))+((ep*(Rog(1))*(W(5,1)))/(dt))+(((Rog(1)*vel(1))/(dz))*(WW(5)))+(((ep*Rog(1)*Dli(5,1))/(0.5*(dz)^2))*(WW(5)));
Aw5(nz,nz-1)=-((((Rog(nz))*(vel(nz)))/(dz))+((ep*(Rog(nz))*(Dli(5,nz)))/((dz)^2)));
Aw5(nz,nz)=((ep*(Rog(nz)))/(dt))+(((RogE)*(vel(nz+1)))/(dz))+((ep*(Rog(nz))*(Dli(5,nz)))/((dz)^2));
bW5(nz)=SW(5,nz)+(((ep*Rog(nz))/(dt))*(W(5,nz)));

for i=2:nz-1
    Aw5(i,i-1)=-((((Rog(i))*(vel(i)))/(dz))+((ep*(Rog(i))*(Dli(5,i)))/((dz)^2)));
    Aw5(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i+1))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(5,i+1)))/((dz)^2))+((ep*(Rog(i))*(Dli(5,i)))/((dz)^2));
    Aw5(i,i+1)=-((ep*(Rog(i+1))*(Dli(5,i+1)))/((dz)^2));
    bW5(i)=SW(5,i)+(((ep*Rog(i))/(dt))*(W(5,i)));
end

%species equation 6
Aw6(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(2))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(6,2)))/((dz)^2))+((ep*(Rog(1))*(Dli(6,1)))/(0.5*(dz)^2));
Aw6(1,2)=-((ep*(Rog(2))*(Dli(6,2)))/((dz)^2));
bW6(1)=(SW(6,1))+((ep*(Rog(1))*(W(6,1)))/(dt))+(((Rog(1)*vel(1))/(dz))*(WW(6)))+(((ep*Rog(1)*Dli(6,1))/(0.5*(dz)^2))*(WW(6)));
Aw6(nz,nz-1)=-((((Rog(nz))*(vel(nz)))/(dz))+((ep*(Rog(nz))*(Dli(6,nz)))/((dz)^2)));
Aw6(nz,nz)=((ep*(Rog(nz)))/(dt))+(((RogE)*(vel(nz+1)))/(dz))+((ep*(Rog(nz))*(Dli(6,nz)))/((dz)^2));
bW6(nz)=SW(6,nz)+(((ep*Rog(nz))/(dt))*(W(6,nz)));

for i=2:nz-1
    Aw6(i,i-1)=-((((Rog(i))*(vel(i)))/(dz))+((ep*(Rog(i))*(Dli(6,i)))/((dz)^2)));
    Aw6(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i+1))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(6,i+1)))/((dz)^2))+((ep*(Rog(i))*(Dli(6,i)))/((dz)^2));
    Aw6(i,i+1)=-((ep*(Rog(i+1))*(Dli(6,i+1)))/((dz)^2));
    bW6(i)=SW(6,i)+(((ep*Rog(i))/(dt))*(W(6,i)));
end

%species equation 7
Aw7(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(2))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(7,2)))/((dz)^2))+((ep*(Rog(1))*(Dli(7,1)))/(0.5*(dz)^2));
Aw7(1,2)=-((ep*(Rog(2))*(Dli(7,2)))/((dz)^2));
bW7(1)=(SW(7,1))+((ep*(Rog(1))*(W(7,1)))/(dt))+(((Rog(1)*vel(1))/(dz))*(WW(7)))+(((ep*Rog(1)*Dli(7,1))/(0.5*(dz)^2))*(WW(7)));
Aw7(nz,nz-1)=-((((Rog(nz))*(vel(nz)))/(dz))+((ep*(Rog(nz))*(Dli(7,nz)))/((dz)^2)));
Aw7(nz,nz)=((ep*(Rog(nz)))/(dt))+(((RogE)*(vel(nz+1)))/(dz))+((ep*(Rog(nz))*(Dli(7,nz)))/((dz)^2));
bW7(nz)=SW(7,nz)+(((ep*Rog(nz))/(dt))*(W(7,nz)));

for i=2:nz-1
    Aw7(i,i-1)=-((((Rog(i))*(vel(i)))/(dz))+((ep*(Rog(i))*(Dli(7,i)))/((dz)^2)));
    Aw7(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i+1))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(7,i+1)))/((dz)^2))+((ep*(Rog(i))*(Dli(7,i)))/((dz)^2));
    Aw7(i,i+1)=-((ep*(Rog(i+1))*(Dli(7,i+1)))/((dz)^2));
    bW7(i)=SW(7,i)+(((ep*Rog(i))/(dt))*(W(7,i)));
end
%Solve the mass fraction global system
W(1,1:nz)=sparse(Aw1)\sparse(bW1');
W(2,1:nz)=sparse(Aw2)\sparse(bW2');
W(4,1:nz)=sparse(Aw4)\sparse(bW4');
W(5,1:nz)=sparse(Aw5)\sparse(bW5');
W(6,1:nz)=sparse(Aw6)\sparse(bW6');
W(7,1:nz)=sparse(Aw7)\sparse(bW7');

%total mass fraction
Wt=W(1,:)+W(2,:)+W(3,:)+W(4,:)+W(5,:)+W(6,:)+W(7,:)
%   Species1PlusBoundary+Species2PlusBoundary+Species3PlusBoundary+Species4PlusBoundary+Species5PlusBoundary+Species6PlusBoundary+Species7PlusBoundary;
%  Conversion from mass to mole fraction
ytb=(W(1,:)./M(1))+(W(2,:)./M(2))+(W(3,:)./M(3))+(W(4,:)./M(4))+(W(5,:)./M(5))+(W(6,:)./M(6))+(W(7,:)./M(7))
% (W(1,:)./M(1))*Wt+(W(2,:)./M(2))*Wt+(W(3,:)./M(3))*Wt+(W(4,:)./M(4))*Wt+(W(5,:)./M(5))+(W(6,:)./M(6))+(W(7,:)./M(7))
ytt=1./ytb;

%mole fraction conversion
y(1,1:nz)=(W(1,1:nz)./(M(1)./ytt))
y(2,1:nz)=(W(2,1:nz)./(M(2)./ytt))
y(3,1:nz)=yw(3);
y(4,1:nz)=(W(4,1:nz)./(M(4)./ytt))
y(5,1:nz)=(W(5,1:nz)./(M(5)./ytt))
y(6,1:nz)=(W(6,1:nz)./(M(6)./ytt))
y(7,1:nz)=(W(7,1:nz)./(M(7)./ytt))

%total mole fraction
yt=y(1,:)+y(2,:)+y(3,:)+y(4,:)+y(5,:)+y(6,:)+y(7,:)

%For the convenience of curve plotting (32,1) each
Species1PlusBoundary=[yw(1) y(1,:) y(1,nz)]';
Species2PlusBoundary=[yw(2) y(2,1:nz) y(2,nz)]';
Species3PlusBoundary=[yw(3) y(3,1:nz) y(3,nz)]';
Species4PlusBoundary=[yw(4) y(4,1:nz) y(4,nz)]';
Species5PlusBoundary=[yw(5) y(5,1:nz) y(5,nz)]';
Species6PlusBoundary=[yw(6) y(6,1:nz) y(6,nz)]';
Species7PlusBoundary=[yw(7) y(7,1:nz) y(7,nz)]';

%Average molecular weight of the mixture in Kg/kmol
for is=1:nz
    Mavg=AvgMolarWt(y(1,:),y(2,:),y(3,:),y(4,:),y(5,:),y(6,:),y(7,:));
end