
tic

Parameters
Params_Init
% PT is failing with the following error
% ?error: 'AvgMolarWt' undefined near line 32 column 10
% error: called from
%    PT at line 32 column 9
%    Anita at line 6 column 1
PT 

toc

figure  (1)
plot(z',T)
xlabel('reactor length (m)') 
ylabel('Pressure (Pa)')
% plot(MeshCell(:,2),MatrixTemp(:,2),'r--o',...
%     MeshCell(:,2),MatrixTemp(:,timestep*(1/4)),'b--o',...
%     MeshCell(:,2),MatrixTemp(:,timestep*(1)),'k--o');
title('Temperature Profile Along Axis'); 
ylabel('Temperature, K'); xlabel('Axial Position, m');
ytickformat('%.1f')