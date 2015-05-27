function [nu,pu]=npcdkp(n,p)

global semiNodes;
global mun mup;
global ni;
global doping;

nu=n;
pu=p;

%munp = mun/mup;
munp = 1;

for i = 1:length(semiNodes)
    nc = n(semiNodes(i));
    pc = p(semiNodes(i));
    if nc<0 && pc<0
     nu(semiNodes(i)) = abs(pc/munp);
     pu(semiNodes(i)) = abs(nc*munp); 
    elseif nc<0 && pc>0
     pu(semiNodes(i)) = pc+abs(nc*munp); 
     nu(semiNodes(i)) = 0*ni^2/pu(semiNodes(i));
    elseif nc>0 && pc<0
     nu(semiNodes(i)) = nc+abs(pc/munp);
     pu(semiNodes(i)) = 0*ni^2/nu(semiNodes(i));
    end
end
