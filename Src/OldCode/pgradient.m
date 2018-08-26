%pgradient is the pressure gradient obtained from Ergun's equation
function Kp= pgradient(mu,uz,ep,dp,Rog)
%     f =((150*mu*uz*((1-ep)^2))/((dp^2)*(ep^3)))+((1.75*Rog*(uz^2)*(1-ep))./dp*(ep^3));
kD =((150*mu*((1-ep)^2))/((dp^2)*(ep^3)))
kV=((1.75*Rog*(1-ep))./dp*(ep^3));
Kp=1./(kD+kV*uz)
end