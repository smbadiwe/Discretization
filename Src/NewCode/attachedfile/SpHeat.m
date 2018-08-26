%Hi is the integral of cpdT
 function [cpg1,Hi]=SpHeat(nz,A,B,C,D,E,T,Tw)
    Hi(1:nz)=0;
    Tao(1:nz)=0;
   Tw;
 Tao=T-Tw;
    cpg1=A+(B.*T)+(C*(T.^2))+(D*(T.^3))+(E*(T.^4));
    Hi=A*Tao+((B/2)*(Tao.^2))+((C/3)*(Tao.^3))+((D/4)*(Tao.^4))+((E/5)*(Tao.^5));
 end
 
