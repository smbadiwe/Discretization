classdef TransportProject < handle
   properties
     appData;
     z;
     dz;
     zEdge;
     nz;

     shcConst; %Constant for specific heat capacity[J/molK] combined.
     Press;
     T;
     T0;
   end
   
   methods
     % Constructor
     function obj = TransportProject()
        appData = AppData;
        nz = appData.nz;
        dz = appData.L/nz;
        % Define a row vector z upfront instead of resizing repeatedly
        z = zeros(1, nz+2); 
        z(1)=0; 
        z(2)=dz/2;  
        z(3:nz+1)=(z(2)+(1:nz-1)*dz);    
        z(nz+2)=z(nz+1)+dz/2;   
        zEdge(1:nz+1) = (0:nz)*dz;
        Press(1:nz) = appData.pE;
        
        obj.zEdge = zEdge;
        obj.z = z;
        obj.nz = nz;
        obj.dz = dz;
        obj.T = appData.T;
        obj.T0 = appData.Tambient;
        obj.Press = Press;
        obj.appData = appData;
        obj.shcConst = [appData.A; appData.B; appData.C; appData.D; appData.E];
        
        % obj.Execute();
     end
    
     function Execute(obj)
        appData = obj.appData;
        nt = appData.nt;
        nz = obj.nz;
        ns = appData.ns;
        Dm = appData.Dm;
        M = appData.M;
        yw = appData.yw;
        
        % Initial values
        %==================
        Mavg = obj.AvgMolarWt();
        y = zeros(nz, ns);
        for is=1:ns
          y(1:nz,is)=yw(is);
        end
        
        W = zeros(nz, ns);
        WW = zeros(ns, 1);
        %Mass fraction of species
        for is=1:ns
          WW(is)=(M(is)*yw(is))./Mavg';
        end
        W(1:nz,1:ns)=WW(ns);
        
        velAvg = appData.uz_in;
        Hrxno = obj.getStandardHeatOfReaction();
        SE = zeros(1,nz); 
        SW = zeros(nz, ns);
        vel = zeros(1, nz+1);
        for it=1:nt
          Rog = obj.getRog(Mavg);
          Kp = obj.GetKp(velAvg, Rog);
          [ApTT, bp] = obj.DensityEquation(Rog, Kp);
          [cpg, H1] = obj.getCpg(Mavg);
          
          vel = obj.getVelocity(Kp);
          [ATT, bT] = obj.EnergyEquation(Rog, vel, cpg, SE);
          
          AG = [ApTT;ATT];
          bG = [bp,bT];
          
          %Solve global system (Solve P and T simultaneously)
          % ISSUE: Sometimes, this 'warns' that matrix is singular (so no inverse)
          PT = sparse(AG)\sparse(bG');
          
          % UPDATES
          obj.Press(1:nz) = PT(1:nz);
          obj.T(1:nz) = PT(nz+1:2*nz);
          
          Dli = obj.GetAxialDispersionCoefficient(velAvg);
          [Aw, bW] = obj.SpeciesEquation(vel, Rog, Dli, SW, W, WW);
          
          %Update the average velocity[m/s]
          velAvg(1:nz) = (vel(1:nz)+vel(2:nz+1))/2; % This seens fishy
          Hrxn = obj.getHeatOfReaction(Hrxno, H1);
          % TODO: Update y
          RR = obj.ReactionRate(y);
          SE = obj.GetSourceTermForTemperature(RR, Hrxn);
          SW = obj.GetSourceTermForSpecies(RR);
        end
        
        % Now show plot
        obj.ShowPlots(vel);
     end
     
     function ShowPlots(obj, vel)
       
        appData = obj.appData;
        nz = obj.nz;
        PressPlusBoundary=[appData.pW obj.Press appData.pE]';
        TempPlusBoundary=[appData.Tw obj.T obj.T(nz)]';  
        figure  (1)
        plot(obj.z',PressPlusBoundary)
        xlabel('reactor length (m)') 
        ylabel('Pressure (Pa)')
        axis([0 0.045 1E5 1.12E5])
        figure (2)
        plot(obj.zEdge',vel)
        xlabel('reactor length (m)') 
        ylabel('velocity (m/s)') 
        axis([0 0.045 0.0 inf])
        figure (3)
        plot(obj.z',TempPlusBoundary)
        xlabel('reactor length (m)') 
        ylabel(' Temperature (K)') 
        axis([0 0.045 772 inf])
     end
     function [Aw, bW] = SpeciesEquation(obj, vel, Rog, Dli, SW, W, WW)
        appData = obj.appData;
        nz = obj.nz;
        dt = appData.dt;
        dz = obj.dz;
        ep = appData.ep;
        
        Aw = zeros(nz);
        bW = zeros(nz, 1);
        % left end
        Aw(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(1))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(2,1)))/((dz)^2))+((ep*(Rog(1))*(Dli(1,1)))/(0.5*(dz)^2));
        Aw(1,2)=-((ep*(Rog(2))*(Dli(2,1)))/((dz)^2));
        bW(1)=(SW(1,1))+(((ep*Rog(1))/(dt))*(W(1,1)))+(((Rog(1)*vel(1))/(dz))*(WW(1)))+(((ep*Rog(1)*Dli(1,1))/(0.5*(dz)^2))*(WW(1)));  
        % right end
        Aw(nz,nz-1)=-((((Rog(nz-1))*(vel(nz)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,1)))/((dz)^2)));
        Aw(nz,nz)=((ep*(Rog(nz)))/(dt))+(((Rog(nz))*(vel(nz+1)))/(dz))+((ep*(Rog(nz-1))*(Dli(nz-1,1)))/((dz)^2));
        bW(nz)=SW(nz,1)+(((ep*Rog(nz))/(dt))*(W(nz,1)));
        % middle
        for i=2:nz-1
          Aw(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(i+1,1)))/((dz)^2))+((ep*(Rog(i-1))*(Dli(i-1,1)))/((dz)^2));
          Aw(i,i-1)=-((((Rog(i-1))*(vel(i)))/(dz))+((ep*(Rog(i-1))*(Dli(i-1,1)))/((dz)^2)));
          Aw(i,i+1)=-((ep*(Rog(i+1))*(Dli(i+1,1)))/((dz)^2));
          bW(i)=SW(i,1)+(((ep*Rog(i))/(dt))*(W(i,1)));
        end
        
        % We're probably not done here and these outpots are probably not what we need
     end
     
     function [cpg, H1] = getCpg(obj, Mavg)
        appData = obj.appData;
        [cP, H1] = obj.SpHeat();

        yw = appData.yw;
        %cpmolmix units KJ/KmolK
        % cPmolmix=(abs(((y(:,7).*cP(:,7))+(y(:,6).*cP(:,6))+(y(:,5).*cP(:,5)))-(y(:,4).*cP(:,4))-(y(:,3).*cP(:,3))-(y(:,2).*cP(:,2))-(y(:,1).*cP(:,1))))
        cPmolmix=(abs(((yw(7).*cP(:,7))+(yw(6).*cP(:,6))+(yw(5).*cP(:,5)))-(yw(4).*cP(:,4))-(yw(3).*cP(:,3))-(yw(2).*cP(:,2))-(yw(1).*cP(:,1))))
        
        cpg=((cPmolmix*1000)./Mavg); %units in J/KgK
     end
     
     function [ATT, bT] = EnergyEquation(obj, Rog, vel, cpg, SE)
        % Energy Equation
        appData = obj.appData;
        ep = appData.ep;
        cpp = appData.cpp;
        Rop = appData.Rop;
        k = appData.k;
        T = appData.T;
        Tw = appData.Tw;
        nz = obj.nz;
        dt = appData.dt;
        dz = obj.dz;
        
        %1st term in temperature equation
        a = (1 - ep) * cpp * Rop;

        AT = zeros(nz);
        bT = zeros(1, nz);
        % left end
        AT(1,1)=((a+(ep*(cpg(1))*(Rog(1))))/(dt))+(((Rog(1))*(cpg(1))*(vel(2)))/(dz))+((3*k)/((dz)^2));
        AT(1,2)=-(k/((dz)^2));
        bT(1)=(SE(1))+(((a+(ep*(cpg(1))*(Rog(1))))/dt)*(T(1)))+((((Rog(1))*(cpg(1))*(vel(1)))/(dz))*(Tw))+((k/(0.5*((dz)^2)))*(Tw));
        % right end
        AT(nz,nz)=((a+(ep*(cpg(nz))*(Rog(nz))))/(dt))+(((Rog(nz))*(cpg(nz))*(vel(nz+1)))/(dz))+(k/((dz)^2));
        AT(nz,nz-1)=-((((Rog(nz-1))*(cpg(nz-1))*(vel(nz)))/(dz))+(k/((dz)^2)));
        bT(nz)=(SE(nz))+(((a+(ep*(Rog(nz))*(cpg(nz))))/(dt))*(T(nz)));
        % in-between
        for iT=2:nz-1
            AT(iT,iT)=((a+(ep*(cpg(iT))*(Rog(iT))))/(dt))+(((Rog(iT))*(cpg(iT))*(vel(iT+1)))/(dz))+((2*k)/((dz)^2));
            AT(iT,iT-1)=-((((Rog(iT-1))*(cpg(iT-1))*(vel(iT)))/(dz))+(k/((dz)^2)));
            AT(iT,iT+1)=-(k/((dz)^2));
            bT(iT)=(SE(iT))+(((a+(ep*(cpg(iT))*(Rog(iT))))/(dt))*(T(iT))) ;
        end
        %AT
        Atl=zeros(nz);
        ATT=[Atl,AT];
     end
     
     function [ApTT, bp] = DensityEquation(obj, Rog, Kp)
        %Density Equation
        appData = obj.appData;
        pW = appData.pW;
        pE = appData.pE;
        ep = appData.ep;
        nz = obj.nz;
        dt = appData.dt;
        dz = obj.dz;
        Ap = zeros(nz);
        bp = zeros(1, nz);
        
        % left end
        Ap(1,1)=((ep*(Rog(1)))/((obj.Press(1))*(dt)))-(((Kp(1))*(Rog(1)))/(0.5*((dz)^2)))-(((Kp(2))*(Rog(2)))/((dz)^2));
        Ap(1,2)=(((Kp(2))*(Rog(1)))/((dz)^2));
        bp(1)=-(((Kp(1))*(Rog(1)))/(0.5*(dz^2)))*(pW);
        % right end
        Ap(nz,nz-1)=(((Kp(nz-1))*(Rog(nz-1)))/((dz)^2));
        Ap(nz,nz)=((ep*(Rog(nz)))/((obj.Press(nz))*(dt)))-(((Kp(nz-1))*(Rog(nz-1)))/((dz)^2))-(((Kp(nz))*(Rog(nz)))/(0.5*(dz)^2));
        bp(nz)=-(((Kp(nz))*(Rog(nz))*pE)/(0.5*((dz)^2)));
        % in-between
        for ip=2:nz-1
            Ap(ip,ip)=((ep*(Rog(ip)))/((obj.Press(ip))*dt))-(((Kp(ip-1))*(Rog(ip-1)))/(dz^2))-(((Kp(ip+1))*(Rog(ip+1)))/((dz)^2));
            Ap(ip,ip+1)=(((Kp(ip+1))*(Rog(ip+1)))/((dz)^2));
            Ap(ip,ip-1)=(((Kp(ip-1))*(Rog(ip-1)))/((dz)^2));
            bp(ip)=0;
        end
        
        Ct(1:nz)=-(ep*(Rog(1:nz)))./(obj.T(1:nz)*dt);
        Ap1 = eye(nz);
        ApT = zeros(nz);
        for i=1:nz
            ApT(i,i)=Ct(i).*Ap1(i,i);
        end
        ApTT=[Ap,ApT];
     end
  
     function [Kp] = GetKp(obj, velAvg, Rog)
        appData = obj.appData;
        nz = obj.nz;
        kD = appData.kD;
        ep = appData.ep;
        dp = appData.dp;

        %kinetic pressure drop [kg/m4]
        kV=(1.75*Rog(1:nz)*(1-ep))./dp*(ep^3); 
        % Kp = 1./(kD(1:nz)+(kV.*(velAvg(1:nz))));
        Kp = 1./(kD+(kV.*(velAvg)));
     end
     
     function [cP, H1] = SpHeat(obj)
       %Hi is the integral of cpdT
        appData = obj.appData;
        ns = appData.ns;
        nz = obj.nz;
        T = obj.T;
        T0 = obj.T0;
        Tao = T-T0;
        
        H1 = zeros(nz, ns); %(J/mol)
        cP = zeros(nz, ns); %units in [J/molK] and [J/KgK]
        for is=1:ns
          % cpg1=A+(B.*T)+(C*(T.^2))+(D*(T.^3))+(E*(T.^4));
          % Hi=A*Tao+((B/2)*(Tao.^2))+((C/3)*(Tao.^3))+((D/4)*(Tao.^4))+((E/5)*(Tao.^5));
          cpg1 = zeros(1, nz);
          Hi = zeros(1, nz);
          A2E = obj.shcConst(:,is);
          for j=1:5
            cpg1 = cpg1 + (A2E(j) * T.^(j - 1));
            Hi = Hi + ((A2E(j)/j) * Tao.^j);
          end
          
          cP(1:nz,is) = cpg1;
          H1(1:nz,is) = Hi;
        end
        
     end
     
     function [Dli] = GetAxialDispersionCoefficient(obj, velAvg)
        %Axial dispersion coefficient in m2/s
        appData = obj.appData;
        nz = obj.nz;
        Dm = appData.Dm;
        dp = appData.dp;
        ns = appData.ns;
        Dli = zeros(nz, ns);
        for is=1:ns
          % Dli(1:nz,is)=0.73.*(Dm(is))+((0.5.*velAvg(1:nz)*dp)./(1+((9.49.*(Dm(is)))./((velAvg(1:nz)).*dp))));
          Dli(1:nz,is)=0.73.*(Dm(is))+((0.5.*velAvg*dp)./(1+((9.49.*(Dm(is)))./((velAvg).*dp))));
        end
     end
     
     function [SE] = GetSourceTermForTemperature(obj, RR, Hrxn)
       
        appData = obj.appData;
        nz = obj.nz;
        ep = appData.ep;
        nk = appData.nk;
        Rop = appData.Rop;
        %Sourceterm for temperature equation - Possible bug here!
        SE(1:nz) =(((1-ep)*(Rop)*((-1*Hrxn)).*nk.*(RR(1:nz))))*1E6; % in [J/m3s]

     end
     function [SW] = GetSourceTermForSpecies(obj, RR)
       
        %Source term for species equation
        appData = obj.appData;
        nz = obj.nz;
        ns = appData.ns;
        ep = appData.ep;
        nk = appData.nk;
        Rop = appData.Rop;
        M = appData.M;
        v = appData.v;
        
        SW = zeros(nz, ns);
        for is=1:ns
          SW(1:nz,is)=(1-ep)*Rop.*(M(is)).*(nk.*(RR(1:nz)).*(v(is)))
        end
     end
     
     function [RR] = ReactionRate(obj, y)
%Rate of reaction[molC/gcat.s] same as [kmolC/kgcat.s]
        % Reaction rate update[molC/gcat.s]
        appData = obj.appData;
        nz = obj.nz;
        ns = appData.ns;
        An = appData.An;
        Apr = appData.Apr;
        Ru = appData.Ru;
        m = appData.m;
        n = appData.n;
        o = appData.o;
        
        %update of partial pressure
        P = zeros(nz, ns);
        for is=1:ns
          P(1:nz,is)=y(1:nz,is).*(obj.Press(1:nz)');
        end
        P
        y
        obj.Press
        m
        n
        o
        Ru
        Apr
        An
        RR(1:nz)=((Apr)*(exp(An./(Ru*(obj.T(1:nz)'))))).*((P(1:nz,1)).^m).*((P(1:nz,4)).^n).*((P(1:nz,2)).^o);
        RR
     end
     
     function [Hrxn] = getHeatOfReaction(obj, Hrxno, H1)
        % %Heat of reaction
        appData = obj.appData;
        nz = obj.nz;
        v = appData.v;
        Hf = appData.Hf;
        H1_mix(1:nz)=(v(7)*(H1(:,7)))+(v(6)*(H1(:,6)))+((v(5)*H1(:,5)))-(v(4)*(H1(:,4)))-(v(3)*(H1(:,3)))-(v(2)*(H1(:,2)))-(v(1)*(H1(:,1)));
        Hrxn(1:nz)=(Hrxno)+((H1_mix(1:nz))/1000); % in KJ/mol
     end
     
     function [Mavg] = AvgMolarWt(obj)
       % M is Molecular weight of species Kg/kmol; a row vector; 
       % yw is a column vector
        appData = obj.appData;
        Mavg = appData.M * appData.yw';
     end
     
     function [Hrxno] = getStandardHeatOfReaction(obj)
        %Standard heat of reaction[KJ/mol]
        appData = obj.appData;
        v = appData.v;
        Hf = appData.Hf;
        Hrxno=((v(5)*Hf(5))+(v(6)*Hf(6))+(v(7)*Hf(7)))-(v(1)*Hf(1))-(v(2)*Hf(2))-(v(3)*Hf(3))-(v(4)*Hf(4));
     end
     
      
     function [vel] = getVelocity(obj, Kp)
        %velocity.This is from Eq 20 in the model (pdf) 
        appData = obj.appData;
        pE = appData.pE;
        pW = appData.pW;
        nz = appData.nz;
        dz = obj.dz;
        Press = obj.Press;
        kD = appData.kD;
        vel = zeros(1, nz+1);
%         vel(1)=(-1*(Press(1)-pW))/ (kD *(0.5*dz));
%         vel(2:nz)=(-1*(Press(2:nz)-Press(1:nz-1)))/(kD*dz);
%         vel(nz+1)=(-1*(pE-Press(nz)))./(kD *(0.5*dz));
        
        % The 'update' part. Figure out what should happen here
        vel(1)=(-Kp(1))*(Press(1)-pW) / (0.5*dz);
        vel(2:nz)=-(Kp(2:nz)).*(Press(2:nz)-Press(1:nz-1))/dz;
        vel(nz+1)=-(Kp(nz))*(pE-Press(nz))/(0.5*dz);
     end
     
     function [Rog] = getRog(obj, Mavg)
        % density of gas mixture in kg/m3. This is from Eq 20 in the model (pdf)
        
        pE = obj.appData.pE;
        nz = obj.appData.nz;
        R = obj.appData.R;
        T = obj.appData.T;
        Rog = zeros(1, nz+1);
        Rog(1:nz)=((obj.Press(1:nz))*Mavg)/(R*(T(1:nz)));
        Rog(nz+1)=(pE*Mavg)/(R*(T(nz)));
     end
     
    
   end
  
end