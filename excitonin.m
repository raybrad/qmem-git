function X = excitonin(V,n,p)

global ehpgr;
global ehpkd;
global ehpkf;
global ehpkr;
global ni;
global isehpgr;
global Nnode Nlink Nsurf;
global semiNodes;

X = zeros(Nnode,1);

if isehpgr == 1
 ehpkdns = calcukd(V);
 for n1 = semiNodes.'
   ehppn = 1/(1*ehpkdns(n1) + ehpkf);
   X(n1) = ehppn*ehpgr;
 end
end
