function chgngrns = chargenetgr(V,n,p,X)

global ehpgr;
global ehpkd;
global ehpkf;
global ehpkr;
global ni;
global isehpgr;
global Nnode Nlink Nsurf;
global semiNodes;

chgngrns = zeros(Nnode,2);

if isehpgr == 1
 ehpkdns = calcukd(V);
 for n1 = semiNodes.'
   if n(n1)>0 && p(n1) >0
      temnp =1;
   else 
      temnp =0;
   end
   chgngrns(n1,1) = ehpkdns(n1)* X(n1) - temnp*ehpkr*(n(n1)*p(n1)-ni*ni);
   chgngrns(n1,2) = ehpgr - (ehpkf + ehpkdns(n1))* X(n1) + temnp*ehpkr*(n(n1)*p(n1)-ni*ni);
 end
end
