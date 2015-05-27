function dchgngrdns = chargenetgr(V,n,p)

global ehpkd;
global ehpkf;
global ehpkr;
global isehpgr;
global ni;
global Nnode Nlink Nsurf;
global semiNodes;

dchgngrdns = zeros(Nnode,4);

if isehpgr==1
 ehpkdns = calcukd(V);
 for n1 = semiNodes.'
   if n(n1)>0 && p(n1) >0
      temnp =1;
   else 
      temnp =0;
   end
   dchgngrdns(n1,1) = -temnp*ehpkr*p(n1);
   dchgngrdns(n1,2) = -temnp*ehpkr*n(n1);
   dchgngrdns(n1,3) =  ehpkdns(n1);
   dchgngrdns(n1,4) = -(ehpkf + ehpkdns(n1));
 end
end
