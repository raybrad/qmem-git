function [currdl,currdr,currdt]=QM_Calc(V,effcd,currtime,dt)
global nLead QMfeedback;
display(['  Start TD QM per step:']);
tic;


emb2qmb(V,effcd,currtime,dt);

! cp -f embound-v embound ./QMwork;  cd ./QMwork;  /home/hy/sandbox/newLodestar_test/lodestar_ExpIn <in.gstd > qm.log; cp -f qmbound* ../;
%% $HOME/lodestar/bin/lodestartd < input  > qm.log; cp -f qmbound* ../;

if (nLead>=1)
load 'qmbound';
end

currdl=0.0;
currdr=0.0;
currdt=0.0;
if (QMfeedback==0)
currdl=0.0;
currdr=0.0;
currdt=0.0;
display(['  in QMfeedback=0.']);
elseif(QMfeedback==1)
display(['  in QMfeedback=1.']);
if (nLead==1)
currdl= qmbound(1);
elseif(nLead==2)
currdl= qmbound(1);
currdr=-qmbound(2);
elseif(nLead==3)
currdl= qmbound(1);
currdr=-qmbound(2);
currdt= qmbound(3);
end
end
%aqmbound=(qmbound(1)-qmbound(2))/2;  %A/m^2


toc;
display(['  End TD QM per step.']);
