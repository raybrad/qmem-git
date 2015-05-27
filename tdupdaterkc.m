function [Vn,nn,pn,An,Hn] = tdupdaterk(V,n,p,A,H,dtVp,dtHp,dt,ntp)

global epdf;
global extfelc;
global dirNodes dcVolDirNodes acVolDirNodes;
global nsteps;

timec = (ntp-1)*dt;
timeh = timec + 0.5*dt;
timee = timec + dt;

switch extfelc
   case 1
     dtVp(dirNodes) = acVolDirNodes*exp(- timec/epdf)/epdf;
   case 2
     dtVp(dirNodes) = pi*acVolDirNodes*sin(2*pi*timec/epdf)/epdf;
   case 3
     dtVp(dirNodes) = 0.0;
   case 4
     dtVp(dirNodes) = -acVolDirNodes*exp((timec-nsteps*dt)/epdf)/epdf;
   otherwise
     error('undefined bias profile');
end

[dtV1,dtn1,dtp1,dtH1,itNr1] = tdcalupdatec(V,n,p,A,H,dtVp,dtHp);


switch extfelc
   case 1
     dtVp(dirNodes) = acVolDirNodes*exp(- timeh/epdf)/epdf;
   case 2
     dtVp(dirNodes) = pi*acVolDirNodes*sin(2*pi*timeh/epdf)/epdf;
   case 3
     dtVp(dirNodes) = 0.0;
   case 4
     dtVp(dirNodes) = -acVolDirNodes*exp((timeh-nsteps*dt)/epdf)/epdf;
   otherwise
     error('undefined bias profile');
end

Vh = V + dtV1*dt/2;
nh = n + dtn1*dt/2;
ph = p + dtp1*dt/2;
Hh = H + dtH1*dt/2;
Ah = A +   H *dt/2;

[dtV2,dtn2,dtp2,dtH2,itNr2] = tdcalupdatec(Vh,nh,ph,Ah,Hh,dtVp,dtHp);

Vh1 = V + dtV2*dt/2;
nh1 = n + dtn2*dt/2;
ph1 = p + dtp2*dt/2;
Hh1 = H + dtH2*dt/2;
Ah1 = A +   Hh*dt/2;

[dtV3,dtn3,dtp3,dtH3,itNr3] = tdcalupdatec(Vh1,nh1,ph1,Ah1,Hh1,dtVp,dtHp);

switch extfelc
   case 1
     dtVp(dirNodes) = acVolDirNodes*exp(- timee/epdf)/epdf;
   case 2
     dtVp(dirNodes) = pi*acVolDirNodes*sin(2*pi*timee/epdf)/epdf;
   case 3
     dtVp(dirNodes) = 0.0;
   case 4
     dtVp(dirNodes) = -acVolDirNodes*exp((timee-nsteps*dt)/epdf)/epdf;
   otherwise
     error('undefined bias profile');
end

Ve = V + dtV3*dt;
ne = n + dtn3*dt;
pe = p + dtp3*dt;
He = H + dtH3*dt;
Ae = A + Hh1 *dt;

[dtV4,dtn4,dtp4,dtH4,itNr4] = tdcalupdatec(Ve,ne,pe,Ae,He,dtVp,dtHp);

dtVt = (dtV1+2*dtV2+2*dtV3+dtV4)/6;
dtnt = (dtn1+2*dtn2+2*dtn3+dtn4)/6;
dtpt = (dtp1+2*dtp2+2*dtp3+dtp4)/6;
dtHt = (dtH1+2*dtH2+2*dtH3+dtH4)/6;
dtAt = (  H +2*Hh  +2*Hh1 +He  )/6;

Vn = V + dt* dtVt;
nn = n + dt* dtnt;
pn = p + dt* dtpt;
Hn = H + dt* dtHt;
An = A + dt* dtAt;


