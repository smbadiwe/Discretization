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
        Press(1:nz) = appData.pW;
        
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
        y = zeros(ns, nz);
        for is=1:ns
          y(is,1:nz)=yw(is);
        end

        W = zeros(ns, nz);
        WW = zeros(1, ns);
        %Mass fraction of species
        for is=1:ns
          WW(is)=(M(is)*yw(is));
          W(is,1:nz)=WW(is);
        end
        
        velAvg = appData.uz_in;
        Hrxno = obj.getStandardHeatOfReaction();
        SE = zeros(1,nz); 
        SW = zeros(ns, nz);
        vel = zeros(1, nz+1);
        SpeciesPlusBoundary = zeros(nz+2, ns);
        RogE = 0;
        for it=1:nt
          Mavg = obj.AvgMolarWt(y);
          [Rog, RogE] = obj.getRog(Mavg, RogE);
          
          Kp = obj.GetKp(velAvg, Rog);
          [ApTT, bp] = obj.DensityEquation(Rog, Kp);
          [cpg, H1] = obj.getCpg(Mavg, y);
          
          vel = obj.getVelocity(Kp);
          [ATT, bT] = obj.EnergyEquation(Rog, vel, cpg, SE);
          
          AG = [ApTT;ATT];
          bG = [bp,bT];
          
          %Solve global system (Solve P and T simultaneously)
          PT = sparse(AG)\sparse(bG');
          
          % UPDATES
          obj.Press(1:nz) = PT(1:nz);
          obj.T(1:nz) = PT(nz+1:2*nz);
          
          [W, y, SpeciesPlusBoundary] = obj.SpeciesEquation(vel, Rog, RogE, SW, W, WW, y, it == nt);
          
          %Update the average velocity[m/s]
          velAvg(1:nz) = (vel(1:nz)+vel(2:nz+1))/2;
          Hrxn = obj.getHeatOfReaction(Hrxno, H1, y);
          
          RR = obj.ReactionRate(y);
          SE = obj.GetSourceTermForTemperature(RR, Hrxn);
          SW = obj.GetSourceTermForSpecies(RR);
          
        end
        
        Conc = obj.GetConcentration();
        % Now show plot
        obj.ShowPlots(vel, SpeciesPlusBoundary, Conc);
     end
     
     function ShowPlots(obj, vel, SpeciesPlusBoundary, Conc)
       
        appData = obj.appData;
        yw = appData.yw;
        nz = obj.nz;
        z = obj.z;
        PressPlusBoundary=[appData.pW obj.Press appData.pE];
        TempPlusBoundary=[appData.Tw obj.T obj.T(nz)];  
        
        %For the convenience of curve plotting
        Species1PlusBoundary=SpeciesPlusBoundary(:, 1);
        sizeBD = size(Species1PlusBoundary)
        Species2PlusBoundary=SpeciesPlusBoundary(:, 2);
        Species3PlusBoundary=SpeciesPlusBoundary(:, 3);
        Species4PlusBoundary=SpeciesPlusBoundary(:, 4);
        Species5PlusBoundary=SpeciesPlusBoundary(:, 5);
        Species6PlusBoundary=SpeciesPlusBoundary(:, 6);
        Species7PlusBoundary=SpeciesPlusBoundary(:, 7);

        figure (1)
        plot(z,PressPlusBoundary)
        xlabel('reactor length (m)') 
        ylabel('Pressure (Pa)')
        axis([0 0.045 1E5 1.5E5])
        figure (2)
        plot(obj.zEdge',vel)
        xlabel('reactor length (m)') 
        ylabel('velocity (m/s)') 
        axis([0 0.045 0.0 inf])
        figure (3)
        plot(obj.z',TempPlusBoundary)
        xlabel('reactor length (m)') 
        ylabel(' Temperature (K)') 
        axis([0 0.045 773 inf])
        figure (4)
        plot(z',Species1PlusBoundary,Species2PlusBoundary,Species3PlusBoundary,Species4PlusBoundary,Species5PlusBoundary,Species6PlusBoundary,Species7PlusBoundary)
        plot(z',Species1PlusBoundary,'',z',Species2PlusBoundary,'',z',Species3PlusBoundary,'',z',Species4PlusBoundary,'',z',Species5PlusBoundary,'',z',Species6PlusBoundary,'',z',Species7PlusBoundary,'')
        xlabel('reactor length (m)') 
        ylabel('mole fraction')
        legend('Species1PlusBoundary','Species2PlusBoundary','Species3PlusBoundary','Species4PlusBoundary','Species5PlusBoundary','Species6PlusBoundary','Species7PlusBoundary')
        axis([0 0.045 0 1])
        figure (5)
        plot(z,Conc)
        xlabel('reactor length (m)') 
        ylabel(' concentration(kmol/m3)') 
        axis([0 0.045 0 inf])

     end
     
     function [Conc] = GetConcentration(obj)
        appData = obj.appData;
        pW = appData.pW;
        T = obj.T;
        Tw = appData.Tw;
        Ru = appData.Ru;
        nz = obj.nz;
        % Concentration update
        ConcW = (pW/(Ru*Tw))/1000;
        Conc1 = (obj.Press(1:nz)./(Ru.*T(1:nz)))/1000;
        ConcE = (pW/(Ru*T(nz)))/1000;
        Conc = [ConcW Conc1 ConcE];
     end
     
     function [W, y, SpeciesPlusBoundary] = SpeciesEquation(obj, vel, Rog, RogE, SW, W, WW, y, lastRun)
        appData = obj.appData;
        nz = obj.nz;
        dt = appData.dt;
        ns = appData.ns;
        dz = obj.dz;
        ep = appData.ep;
        M = appData.M;
        yw = appData.yw;
        Dm = appData.Dm;
        dp = appData.dp;
        
        % Wt = zeros(1, nz);
        ytb = zeros(1, nz);
        Dli = zeros(ns, nz+1);
        %Axial dispersion coefficient in m2/s
        for is=1:ns
            Dli(is,1:nz+1)=(0.73.*(Dm(is)))+((0.5*(vel(1:nz+1))*dp)./(1+(9.49.*(Dm(is)))./((vel(1:nz+1)).*dp)));
        end
        
        for nsIndex=1:ns % from nsIndex = 1 to 7
            if nsIndex == 3
                W(nsIndex,1:nz) = WW(nsIndex);
            else
                Aw = zeros(nz);
                bW = zeros(1, nz);
                % left end====See as attached to CalcSpecies
                Aw(1,1)=((ep*(Rog(1)))/(dt))+(((Rog(2))*(vel(2)))/(dz))+((ep*(Rog(2))*(Dli(nsIndex,2)))/((dz)^2))+((ep*(Rog(1))*(Dli(nsIndex,1)))/(0.5*(dz)^2));
                Aw(1,2)=-((ep*(Rog(2))*(Dli(nsIndex,2)))/((dz)^2));
                bW(1)=(SW(nsIndex,1))+((ep*(Rog(1))*(W(nsIndex,1)))/(dt))+(((Rog(1)*vel(1))/(dz))*(WW(nsIndex)))+(((ep*Rog(1)*Dli(nsIndex,1))/(0.5*(dz)^2))*(WW(nsIndex)));

                % right end
                Aw(nz,nz-1)=-((((Rog(nz))*(vel(nz)))/(dz))+((ep*(Rog(nz))*(Dli(nsIndex,nz)))/((dz)^2)));
                Aw(nz,nz)=((ep*(Rog(nz)))/(dt))+(((RogE)*(vel(nz+1)))/(dz))+((ep*(Rog(nz))*(Dli(nsIndex,nz)))/((dz)^2));
                bW(nz)=SW(nsIndex,nz)+(((ep*Rog(nz))/(dt))*(W(nsIndex,nz)));

                % middle
                for i=2:nz-1
                    Aw(i,i-1)=-((((Rog(i))*(vel(i)))/(dz))+((ep*(Rog(i))*(Dli(nsIndex,i)))/((dz)^2)));
                    Aw(i,i)=((ep*(Rog(i)))/(dt))+(((Rog(i+1))*(vel(i+1)))/(dz))+((ep*(Rog(i+1))*(Dli(nsIndex,i+1)))/((dz)^2))+((ep*(Rog(i))*(Dli(nsIndex,i)))/((dz)^2));
                    Aw(i,i+1)=-((ep*(Rog(i+1))*(Dli(nsIndex,i+1)))/((dz)^2));
                    bW(i)=SW(nsIndex,i)+(((ep*Rog(i))/(dt))*(W(nsIndex,i)));
                end

                %Solve the mass fraction global system
                W(nsIndex,1:nz) = sparse(Aw)\sparse(bW');

            end
            
            %summing for total mass fraction - not used, so not computed
            % Wt = Wt + W(nsIndex,:);
            
            %Conversion from mass to mole fraction 
            ytb = ytb + (W(nsIndex,:)./M(nsIndex));
        end
        
        ytt = 1./ytb;
        
        for nsIndex=1:ns 
            if nsIndex == 3
                y(nsIndex,1:nz) = yw(nsIndex);
            else
                y(nsIndex,1:nz) = (W(nsIndex,1:nz)./(M(nsIndex)./ytt));
            end
        end
        SpeciesPlusBoundary = [];
        if lastRun % we only need to compute SpeciesPlusBoundary in the last run for graph plotting
            SpeciesPlusBoundary = zeros(nz+2, ns);
            for nsIndex=1:ns 
                SpeciesPlusBoundary(:, nsIndex) = [yw(nsIndex) y(nsIndex,1:nz) y(nsIndex,nz)]';
            end
        end
     end
     
     function [cpg, H1] = getCpg(obj, Mavg, y)
        [cP, H1] = obj.SpHeat();

        %cpmolmix units KJ/KmolK
        cPmolmix=(abs(((y(7,:).*cP(7,:))+(y(6,:).*cP(6,:))+(y(5,:).*cP(5,:)))-(y(4,:).*cP(4,:))-(y(3,:).*cP(3,:))-(y(2,:).*cP(2,:))-(y(1,:).*cP(1,:))));
        cpg=((cPmolmix*1000)./Mavg); %units in J/KgK
     end
     
     function [ATT, bT] = EnergyEquation(obj, Rog, vel, cpg, SE)
        % Energy Equation
        appData = obj.appData;
        ep = appData.ep;
        cpp = appData.cpp;
        Rop = appData.Rop;
        k1 = appData.k1;
        k2 = appData.k2;
        h = appData.h;
        d = appData.d;
        T = obj.T;
        Tw = appData.Tw;
        nz = obj.nz;
        dt = appData.dt;
        dz = obj.dz;
        
        %1st term in temperature equation
        a = (1 - ep) * cpp * Rop;

        AT = zeros(nz);
        bT = zeros(1, nz);
        % left end
        AT(1,1)=((a+(ep*(cpg(1))*(Rog(1))))/(dt))+(((Rog(2))*(cpg(2))*(vel(2)))/(dz))+((k1)/(0.5*(dz)^2))+((k2)/((dz)^2));
        AT(1,2)=-(k2/((dz)^2));
        bT(1)=(SE(1))+(((4*h)/d)*(Tw-T(1)))+(((a+(ep*(cpg(1))*(Rog(1))))/dt)*(T(1)))+((((Rog(1))*(cpg(1))*(vel(1)))/(dz))*(Tw))+((k1*(Tw))/(0.5*((dz)^2)));
        % right end
        AT(nz,nz)=((a+(ep*(cpg(nz))*(Rog(nz))))/(dt))+(((Rog(nz))*(cpg(nz))*(vel(nz)))/(dz))+(k1/((dz)^2));
        AT(nz,nz-1)=-(((Rog(nz))*(cpg(nz))*(vel(nz)))/(dz))-(k1/((dz)^2));
        bT(nz)=(SE(nz))+(((a+(ep*(Rog(nz))*(cpg(nz))))/(dt))*(T(nz)))+(((4*h)/d)*(Tw-T(nz)));
    
        % in-between
        for iT=2:nz-1
            AT(iT,iT-1)=-(((Rog(iT))*(cpg(iT))*(vel(iT)))/(dz))-(k1/((dz)^2));
            AT(iT,iT)=((a+(ep*(cpg(iT))*(Rog(iT))))/(dt))+(((Rog(iT))*(cpg(iT))*(vel(iT)))/(dz))+((k1)/((dz)^2))+((k2)/((dz)^2));
            AT(iT,iT+1)=-(k2/((dz)^2));
            bT(iT)=(SE(iT))+(((4*h)/d)*(Tw-T(iT)))+(((a+(ep*(cpg(iT))*(Rog(iT))))/(dt))*(T(iT)));
        end
        
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
        Ap(1,1)=((ep*(Rog(1)))/((obj.Press(1))*(dt)))+(((Kp(1))*(Rog(1)))/(0.5*((dz)^2)))+(((Kp(2))*(Rog(2)))/((dz)^2));
        Ap(1,2)=-(((Kp(2))*(Rog(2)))/((dz)^2));
        bp(1)=(((Kp(1))*(Rog(1)))/(0.5*(dz^2)))*(pW);
        % right end
        Ap(nz,nz-1)=-(((Kp(nz-1))*(Rog(nz-1)))/((dz)^2));
        Ap(nz,nz)=((ep*(Rog(nz)))/((obj.Press(nz))*(dt)))+(((Kp(nz-1))*(Rog(nz-1)))/((dz)^2))+(((Kp(nz))*(Rog(nz)))/(0.5*(dz)^2));
        bp(nz)=(((Kp(nz))*(Rog(nz))*pE)/(0.5*((dz)^2)));
        % in-between
        for ip=2:nz-1
            Ap(ip,ip)=((ep*(Rog(ip)))/((obj.Press(ip))*dt))+(((Kp(ip-1))*(Rog(ip-1)))/(dz^2))+(((Kp(ip+1))*(Rog(ip+1)))/((dz)^2));
            Ap(ip,ip+1)=-(((Kp(ip+1))*(Rog(ip+1)))/((dz)^2));
            Ap(ip,ip-1)=-(((Kp(ip-1))*(Rog(ip-1)))/((dz)^2));
            bp(ip)=0;
        end
        
        Ct = -(ep*(Rog(1:nz)))./(obj.T(1:nz)*dt);
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
        
        H1 = zeros(ns, nz); %(J/mol)
        cP = zeros(nz, nz); %units in [J/molK] and [J/KgK]
        for is=1:ns
          % cpg1=A+(B.*T)+(C*(T.^2))+(D*(T.^3))+(E*(T.^4));
          % Hi=A*Tao+((B/2)*(Tao.^2))+((C/3)*(Tao.^3))+((D/4)*(Tao.^4))+((E/5)*(Tao.^5));
          cpg1 = zeros(1, nz);
          Hi = zeros(1, nz);
          A2E = obj.shcConst(:, is);
          for j=1:5
            cpg1 = cpg1 + (A2E(j) * T.^(j - 1));
            Hi = Hi + ((A2E(j)/j) * Tao.^j);
          end
          
          cP(is,1:nz) = cpg1;
          H1(is,1:nz) = Hi;
        end
        
     end
     
     function [SE] = GetSourceTermForTemperature(obj, RR, Hrxn)
       
        appData = obj.appData;
        nz = obj.nz;
        ep = appData.ep;
        nk = appData.nk;
        Rop = appData.Rop;
        % %Sourceterm for temperature equation
        SE(1:nz)=1E6;%((1-ep).*(Rop).*(((-1*Hrxn)).*nk.*(RR(1:nz))));
 
     end
     
     function [SW] = GetSourceTermForSpecies(obj, RR)
       
        %Source term for species equation
        appData = obj.appData;
        nz = obj.nz;
        ns = appData.ns;
        ep = appData.ep;
        % nk = appData.nk;
        Rop = appData.Rop;
        M = appData.M;
        v = appData.v;
        
        SW = zeros(ns, nz);
        multiplier = -1;
        for is=1:ns
            if (is > 4)
                multiplier = 1;
            end
          SW(is,1:nz)=(1-ep)*Rop.*(M(is)).*(RR(1:nz).*(multiplier*v(is))) ;
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
        T = appData.T;
        
        %update of partial pressure
        P = zeros(ns, nz);
        for is=1:ns
          P(is,:)=y(is,:).*(obj.Press(:)');
        end

        %Rate of reaction[molC/gcat.s] same as [kmolC/kgcat.s]
        RR(1:nz)=((Apr*(exp(An./(Ru*(T(1:nz)))))).*((P(1,:)).^m).*((P(4,:)).^n).*((P(2,:)).^o));
     end
     
     function [Hrxn] = getHeatOfReaction(obj, Hrxno, H1, y)
        % %Heat of reaction
        appData = obj.appData;
        nz = obj.nz;
        v = appData.v;
        Hf = appData.Hf;
        H1_mix=((y(7,:).*H1(7,:))+(y(6,:).*H1(6,:))+(y(5,:).*H1(5,:)))-(y(4,:).*H1(4,:))-(y(3,:).*H1(3,:))-(y(2,:).*H1(2,:))-(y(1,:).*H1(1,:));
        %Heat of reaction term
        Hrxn=((Hrxno*1E3)+(2315.98)+((H1_mix)))*1E3;% in J/kmol +(2315.98/1000)
     end
     
     function [Mavg] = AvgMolarWt(obj, y)
       % M is Molecular weight of species Kg/kmol; a row vector; 
       % yw is a column vector
       
       % Mavg=AvgMolarWt(y(1,:),y(2,:),y(3,:),y(4,:),y(5,:),y(6,:),y(7,:));
        appData = obj.appData;
        Mavg = zeros(1, obj.nz);
        for i=1:appData.ns
            Mavg = Mavg + (appData.M(i) * y(i,:));
        end
        % Mavg = appData.M * appData.yw;
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
        
        % The 'update' part. Figure out what should happen here
        vel(1)=(-Kp(1))*(Press(1)-pW) / (0.5*dz);
        vel(2:nz)=-(Kp(2:nz)).*(Press(2:nz)-Press(1:nz-1))/dz;
        vel(nz+1)=-(Kp(nz))*(pE-Press(nz))/(0.5*dz);
     end
     
     function [Rog, RogE] = getRog(obj, Mavg, RogE)
        % density of gas mixture in kg/m3. This is from Eq 20 in the model (pdf)
        pE = obj.appData.pE;
        nz = obj.appData.nz;
        R = obj.appData.R;
        
        %Density Update
        Rog = ((obj.Press(1:nz)).*Mavg)./(R*(obj.T(1:nz)));
        
        if RogE == 0
            RogE = (pE.*Mavg(nz))/(R.*(obj.T(nz)));%same as Rog at nz+1
        end
     end
     
    
   end
  
end