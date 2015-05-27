function X = excitonbc(V,n,p,X)

global ehpgr;
global ehpkd;
global ehpkf;
global ehpkr;
global ni;
global isehpgr;
global Nnode Nlink Nsurf;
global semiNodes;
global bndNodes edgeNodes dirNodes;
global dirSemiNodes isDirSemiNodes;

excitNodes1 = unique([dirSemiNodes;dirNodes]);
excitNodes2 = intersect(semiNodes,excitNodes1);

if isehpgr == 1
 ehpkdns = calcukd(V);
 for n1 = excitNodes2.'
   ehppn = 1/(1*ehpkdns(n1) + ehpkf);
   X(n1) = ehppn*ehpgr;
 end
end
