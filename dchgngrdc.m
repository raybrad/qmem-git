function dchgngrdns = chargenetgr(V,n,p)

global ehpkd;
global ehpkf;
global ehpkr;
global isehpgr;
global ni;
global Nnode Nlink Nsurf;
global semiNodes;

dchgngrdns = zeros(Nnode,2);

if isehpgr==1
 ehpkdns = calcukdc(V);
 ehpp = ehpkdns./(ehpkdns + ehpkf);
 ehpgr2 = zeros(Nnode,1);
 for n1 = semiNodes.'
  if n(n1)>0 && p(n1) >0
   ehpgr2(n1) = -(1-ehpp(n1))*ehpkr;
  end
 end

 dchgngrdns(:,1) = ehpgr2.*p;
 dchgngrdns(:,2) = ehpgr2.*n;
end
