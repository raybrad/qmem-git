function [currdl,currdr,currdt]=QM_Calc(V,effcd,currtime,dt)

display(['  Start TD QM per step:']);
tic;


emb2qmb(V,effcd,currtime,dt);

! cp -f embound-v embound ./QMwork;  cd ./QMwork;  /home/hy/sandbox/newLodestar_test/lodestar_ExpIn <in.gstd > qm.log; cp -f qmbound* ../;
%% $HOME/lodestar/bin/lodestartd < input  > qm.log; cp -f qmbound* ../;

load 'qmbound';

currdl= qmbound(1);
currdr=-qmbound(2);
currdt= qmbound(3);
%aqmbound=(qmbound(1)-qmbound(2))/2;  %A/m^2


toc;
display(['  End TD QM per step.']);
