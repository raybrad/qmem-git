function chgngrns = chargenetgr(V,n,p)

global ehpgr;
global ehpkd;
global ehpkf;
global ehpkr;
global ni;
global isehpgr;
global Nnode Nlink Nsurf;
global semiNodes;

chgngrns = zeros(Nnode,1);

if isehpgr == 1
 ehpkdns = calcukdc(V);
 ehpp = ehpkdns./(ehpkdns + ehpkf);
 ehpgr1 = ehpp*ehpgr;
 ehpgr2 = zeros(Nnode,1);
 for n1 = semiNodes.'
  if n(n1)>0 && p(n1) >0
   ehpgr2(n1) = ehpkr*(1-ehpp(n1))*(n(n1)*p(n1)-ni*ni);
  end
 end
 chgngrns = ehpgr1 - ehpgr2;
end
