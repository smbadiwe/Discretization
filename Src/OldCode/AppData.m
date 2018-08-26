classdef AppData < handle
  properties (SetAccess = private)
      L=0.045; %0.25 % in m
      
      % cpp: number of catalyst pellet[J/KgK]
      cpp=379.64; 
              
      %Constant for specific heat capacity[J/molK]
      A=[9.656      29.526    29.342      33.933    25.399      34.942     27.437];
      B=[4.283E-01 -8.99E-03 -3.54E-03   -8.42E-03  2.0178E-02 -3.996E-02  4.23E-02];
      C=[-2.68E-04  3.81E-05  1.0076E-05  2.99E-05 -3.85E-05    1.92E-04  -1.96E-05];
      D=[3.18E-08  -3.26E-08 -4.31E-09   -1.78E-08  3.19E-08   -1.53E-07   3.997E-09];
      E=[2.77E-11   8.86E-12  2.59E-13    3.69E-12 -8.76E-12    3.93E-11  -2.99E-13];

      % stoichiometric coefficient[-]
      v=[1 0.65 0 0.7 3.2 0.5 2];

      %standard heat of formation of species[kJ/mol]
      Hf=[-582.8 0 0 -241.8 0 -74.85 -393.51];


      %Molecular weight of species Kg/kmol
      M=[69 32 28 18 2 16 44];

      %Diffusion coefficient
      Dm = [3.8E-5  1.03E-04  0  1.14E-04   3.64E-04  1.15E-04  7.94E-05];

      % =================================
      
      %universal gas constant in J/kmolK
      R=8.314*1000; 
      
      Ru=8.314; %[J/molK]
      
      %inlet temeperature in K
      Tambient=298;
      
      Tw=773; %in K
      

      %density of pellet in kg/m3
      Rop=53.5;%8.77;%53.5;
      %effective thermal conductivity in W/mK where k(1) is at the edge
      k=0.062;
      %heat transfer coefficent in W/m2K
      h=14;

      %Boundary & Initial Conditions For Flow and Heat
      pW=1.1E5; %[in Pa]
      pE=1.01325E5; %[Pa}
      
      %Time
      tt=15; %TotalTime in s
      nt=10000; % from 40, we begin to get NaNs in complex number form

      %effectiveness factor for the rate
      nk=0.98;
      
      %%Activation Energy[J/mol] this is because R is also in Kmol
      An=-93540;
      
      % These are set in the constructor
      yw;
      T;
      nz; 
      
      FAO; %inlet flowrate
      ns; % ns: The number of species
      d; % diameter of reactor in m
      dt; % delta time
      m; % Reaction order wrt crude glycerol[-]
      n; % Reaction order wrt water[-]
      o; % Reaction order wrt oxygen[-]
      ep; % porosity
      dp; % diameter of pellets in m
      Apr; % frequency factor  units are [kmolC/kgcat.s.Pa^(m+n+o)
      Area; % Area
      uz_in;
      kD; % viscous pressure drop Units [kg/m4s]
      uv; % viscosity in Pa.s %conversion to Pa.s
        
  end
  
  methods
     % Constructor
     function obj = AppData()
        nz=15; 
        % ns: The number of species
        ns=7;
        %diameter of reactor in m
        d=0.0127;
        dt=1E-04;
        %Reaction order wrt crude glycerol[-]
        m=1.06;
        %Reaction order wrt water[-]
        n=0.56;
        %Reaction order wrt oxygen[-]
        o=1.86;
        % diameter of pellets in m
        dp=0.0008; 
        %inlet flowrate
        FAO=2.10049E-06;

        %Boundary & Initial Conditions For Species
        %Mole fraction of species[mol/mol]
        yw = zeros(1, ns) % ns = 7. Change thi next 2 lines if that changes
        yw(1) = 0.1;  yw(2) = 0.06; yw(3) = 0.23; yw(5) = 0.00001; yw(6) = 0.00001; yw(7) = 0.00001;
        yw(4) = 1-(yw(1)+yw(2)+yw(3)+yw(5)+yw(6)+yw(7));
      
        T(1:nz)=773;
        
        %porosity
        ep=0.38+0.073*(1+((((d/dp)-2)^2)/((d/dp)^2)));

        %frequency factor  units are [kmolC/kgcat.s.Pa^(m+n+o)
        Apr=(((2.135*10^11))/(60*((101325)^(m+n+o))));

        %Area
        Area=(3.142*(d^2))/4;
        uz_in=FAO/Area;
        
        %viscosity in Pa.S
        uv = 0.204;%Pa.s %conversion to Pa.s
        
        %viscous pressure drop Units [kg/m4s]
        kD = (150 * uv * ((1-ep)^2))/((dp^2)*(ep^3)); %229188.56 

        obj.uz_in = uz_in;
        obj.kD = kD;
        obj.Apr = Apr;
        obj.uv = uv;
        obj.ep = ep;
        obj.Area = Area;
        obj.d = d;
        obj.dt = dt;
        obj.ns = ns;
        obj.nz = nz;
        obj.T = T;
        obj.yw = yw;
        obj.dp = dp;
        obj.FAO = FAO;
        obj.m = m;
        obj.n = n;
        obj.o = o;
     end
  end
end
