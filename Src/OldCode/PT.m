R = 8.314; Rop=8.77; cpp=379.64; ns=7; uz_in=0.0663
d=0.0127; dp=0.0008; L=0.045; t=15; nt=1000;dt=t/nt; nz=15;dz=L/nz; pE=101325; Press(1:nz)=pE; uv=0.204;T0=298;T(1:nz)=773; Tao=T-T0;
ep = 0.38 + 0.073 *(1+(((d/dp-2)^2)/((d/dp)^2))); 
ContinuityMatrixLHS = zeros(nz,2*nz); ContinuityMatrixRHS =zeros(nz,1); TempMatrixLHS =zeros(nz,2*nz); TempMatrixRHS =zeros(nz,1);
PressMatrix(1:nz,1:nt)=0; TempMatrix(1:nz,1:nt)=0; MassfracMatrixCG(1:nz,1:nt)=0; MassfracMatrixO2(1:nz,1:nt)=0; MassfracMatrixN2(1:nz,1:nt)=0;MassfracMatrixH2O(1:nz,1:nt)=0;
MassfracMatrixH2(1:nz,1:nt)=0; MassfracMatrixCH4(1:nz,1:nt)=0; MassfracMatrixCO2(1:nz,1:nt)=0;
MatrixVel(1:nz+1,1:nt)=0; MatrixDens(1:nz,1:nt)=0; Matrixcpg(1:nz,1:nt)=0; Matrixkt(1:nz,1:nt)=0; Kp(1:nz)=0; uz(1:nz+1)=0;
 a=(1-ep)*cpp*Rop
for k = 2:nt
for i = 1:nz
  P = PressMatrix(i,k+1);   T = TempMatrix(i,k+1);  yCG=MassfracMatrixCG(i,k+1);  yO2=MassfracMatrixO2(i,k+1);
 yH2O=MassfracMatrixH2O(i,k+1);  yN2=MassfracMatrixN2(i,k+1); yCH4=MassfracMatrixCH4(i,k+1);
 yH2=MassfracMatrixH2(i,k+1);  yCO2=MassfracMatrixCO2(i,k+1); 
 
 im1 = i-1; ip1 = i+1;
 
 uz_in=uz(1);
 
MCG=69; MO2=32; MN2=28; MH2O=18;  MH2=2;  MCH4=16;  MCO2=44;%Molecular weight of species Kg/kmol
%     Rog=MixDensity(Press,Mavg,R,T)
%     cpg = cpmix(nz,yCG,yO2,yH2O,yN2,yCH4,yH2,yCO2,MCG,MO2,MN2,MH2O,MCO2,MH2,MCH4,Mavg)
A(1:ns)=[9.656      29.526    29.342      33.933    25.399      34.942     27.437];
B(1:ns)=[4.283E-01 -8.99E-03 -3.54E-03   -8.42E-03  2.0178E-02 -3.996E-02  4.23E-02];
C(1:ns)=[-2.68E-04  3.81E-05  1.0076E-05  2.99E-05 -3.85E-05    1.92E-04  -1.96E-05];
D(1:ns)=[3.18E-08  -3.26E-08 -4.31E-09   -1.78E-08  3.19E-08   -1.53E-07   3.997E-09];
E(1:ns)=[2.77E-11   8.86E-12  2.59E-13    3.69E-12 -8.76E-12    3.93E-11  -2.99E-13];

 for is=1:ns
 Hi=A(is)*Tao+((B(is)/2).*(Tao.^2))+((C(is)/3).*(Tao.^3))+((D(is)/4).*(Tao.^4))+((E(is)/5)*(Tao.^5));
 end
   Hrxn= Heatrxncalc(T,Hi)
    Mavg=AvgMolarWt(yCG,yO2,yN2,yH2,yH2O,yCO2,yCH4)
    
 kD=(150*uv.*((1-ep)^2))/((dp^2)*(ep^3)); kV=(1.75*Rog*(1-ep))./dp*(ep^3);

uz=Velmix(Press,uz,Kp,kD,kV,nz)

Kp(1:nz)=1./(kD+(kV.*((uz(i)+uz(i+1))/2)))

 uz(1:nz) = MatrixVel(i+1,k-1); Rog(1:nz)= MatrixDens(i,k-1); cpg(1:nz)=Matrixcpg(i,k-1); kt=Matrixkt(i,k-1);
  
SE(1:nz)=0;
 if(i == 1)
     ContinuityMatrixLHS(i,i)=((ep*(Rog(i)))/(Press(i)*(dt)))+((Kp(i)*Rog(i))/(0.5*((dz)^2)))+(((Kp(i+1))*(Rog(i+1)))/((dz)^2))
     ContinuityMatrixLHS(i,i+l)=-(((Kp(i+1))*(Rog(i+1)))/((dz)^2))
     ContinuityMatrixLHS(i,i+nz) =-(ep*(Rog(i))/(T(i)*dt));%(-rho/T*dz/dt)
     ContinuityMatrixRHS(i,1)=((((Kp(i))*(Rog(i)))/(0.5*(dz^2)))*(pW))%((ep*(Rog1))/(dt))
    
     TempMatrixLHS(i,i+nz)=((a+(ep*(cpg(i))*(Rog(i))))/(dt))+(((Rog(i))*(cpg(i))*(uz(2)))/(dz))+((3*kt)/((dz)^2))
     TempMatrixLHS(i,ipl+nz)=-(kt/(dz^2))
     TempMatrixRHS(i,1)=(SE(i))+(((a+(ep*(cpg(i))*(Rog(i))))/dt)*(T(i)))+((((Rog(i))*(cpg(i))*(uz(i)))/(dz))*(Tw))+((kt/(0.5*((dz)^2)))*(Tw))
%       (rho*Cp)*(dz/dt)*T+(-(4*U/d_t)*(T-T_surr))*(dz)-(-(rho1*Cp1*uz1)-cond1/(dz/2))*T_in;
    %Rog is the boundary value;Rog1 is the density value in the 1st
    %cell;Rog2 is the value in the last cell
      elseif (i==nz)
     ContinuityMatrixLHS(i,i-1)=-(((Kp(i-1))*(Rog(i-1)))/((dz)^2))
     ContinuityMatrixLHS(i,i)=((ep*(Rog(i)))/((Press(i))*(dt)))+(((Kp(i-1))*(Rog(i-1)))/((dz)^2))+(((Kp(i))*(Rog(i)))/(0.5*(dz)^2))
     ContinuityMatrixLHS(i,i+nz) =-(ep*(Rog(i)))/(T(i)*dt);%(-rho/T*dz/dt)
     ContinuityMatrixRHS(i,1)=(((Kp(i))*(Rog(i))*pE)/(0.5*((dz)^2)))
    
     TempMatrixLHS(i,i-1+nz)=-((((Rog(i-1))*(cpg(i-1))*(uz(i)))/(dz))+(kt/((dz)^2)))
     TempMatrixLHS(i,i+nz)=((a+(ep*(cpg(i))*(Rog(i))))/(dt))+(((Rog(i))*(cpg(i))*(vel(i+1)))/(dz))+(kt/((dz)^2))
     TempMatrixRHS(i,1)=(SE(i))+(((a+(ep*(Rog(i))*(cpg(i))))/(dt))*(T(i)))
 else
    ContinuityMatrixLHS(i,i-1)=-(((Kp(i-1))*(Rog(i-1)))/((dz)^2))
     ContinuityMatrixLHS(i,i)=((ep*(Rog(i)))/((Press(i))*dt))+(((Kp(i-1))*(Rog(i-1)))/(dz^2))+(((Kp(i))*(Rog(i)))/((dz)^2))
     ContinuityMatrixLHS(i,ipl) =-(((Kp(i))*(Rog(i)))/((dz)^2))
     ContinuityMatrixLHS(i,i+nz) =-(ep*(Rog(i)))/(T(i)*dt)
     ContinuityMatrixRHS(i,1)=0;
    
     TempMatrixLHS(i,iml+nz)=-((((Rog(i-1))*(cpg(i-1))*(uz(i)))/(dz))+(kt/((dz)^2)))
     TempMatrixLHS(i,i+nz)=((a+(ep*(cpg(i))*(Rog(i))))/(dt))+(((Rog(i))*(cpg(i))*(vel(i+1)))/(dz))+((2*kt)/((dz)^2))
      TempMatrixLHS(i,ipl+nz)=-(kt/((dz)^2))
     TempMatrixRHS(i,1)=(SE(i))+(((a+(ep*(cpg(i))*(Rog(i))))/(dt))*(T(i))) 
 end
 
end
LHSMatrixConcatenation = cat(1,ContinuityMatrixLHS,TempMatrixLHS);
RHSMatrixConcatenation = cat(1,ContinuityMatrixRHS,TempMatrixRHS);
PressureTemp = sparse(LHSMatrixConcatenation\RHSMatrixConcatenation);

Press(:,k+1)=PressureTemp(1:nz,1);
T(:,k+1)=PressureTemp(nz+1:2*nz,1);
end