function doping = dopeUniform(semiNodes,QMnodes)

global Nnode Ndp;

doping = zeros(Nnode,1);
%doping(semiNodes) = Ndp * (1+ 5.5 * randn(length(semiNodes),1));
%doping(QMnodes)   = 0;
%doping(semiNodes) = Ndp;
doping(semiNodes) = 0;
%doping(QMnodes)   = 0;
