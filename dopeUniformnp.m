function doping = dopeUniform(semiNodesn,semiNodesp,semiintf,doprf)

global Nnode Ndp Nap ni;

doping = zeros(Nnode,1);
%doping(semiNodes) = Ndp * (1+ 0.05 * randn(length(semiNodes),1));
%doping(QMnodes)   = 0;
doping(semiNodesn) = Ndp - Nap*doprf(semiNodesn);
doping(semiNodesp) = Ndp*doprf(semiNodesp) - Nap;
doping(semiintf) = Ndp - Nap;
%doping(QMnodes)   = 0;
