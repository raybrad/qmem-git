function [Vn,nn,pn,An,Hn] = tdupdaterk(V,n,p,A,H,dtVp,dtHp,dt,ntp)

global epdf;
global extfelc;
global dirNodes dcVolDirNodes acVolDirNodes;
global nsteps;

timec = (ntp-1)*dt;
timeh = timec + 0.5*dt;
timee = timec + dt;

updateTol = 1e-10;
maxNewtonIt = 20;
normUpdate = 1;
itNr = 0;

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

%[dtV1,dtn1,dtp1,dtH1,itNr1] = tdcalupdatec(V,n,p,A,H,dtVp,dtHp);
[dtV1,dtn1,dtp1,dtH1] = tdcalupdates(V,n,p,A,H,dtVp,dtHp);

dtV2 = dtV1;
dtn2 = dtn1;
dtp2 = dtp1;
dtH2 = dtH1;

dtV3 = dtV1;
dtn3 = dtn1;
dtp3 = dtp1;
dtH3 = dtH1;

Hh   = H;
He   = H;

while itNr < maxNewtonIt && normUpdate > updateTol

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

 Vh = V + dt*( 5/24 * dtV1 + 1/3 * dtV2 - 1/24 * dtV3);
 nh = n + dt*( 5/24 * dtn1 + 1/3 * dtn2 - 1/24 * dtn3);
 ph = p + dt*( 5/24 * dtp1 + 1/3 * dtp2 - 1/24 * dtp3);
 Hh = H + dt*( 5/24 * dtH1 + 1/3 * dtH2 - 1/24 * dtH3);
 Ah = A + dt*( 5/24 *    H + 1/3 *   Hh - 1/24 *   He);

% [dtV2t,dtn2t,dtp2t,dtH2t,itNr2] = tdcalupdatec(Vh,nh,ph,Ah,Hh,dtVp,dtHp);
 [dtV2t,dtn2t,dtp2t,dtH2t] = tdcalupdates(Vh,nh,ph,Ah,Hh,dtVp,dtHp);

 rhs_F1  = dtV2t - dtV2;
 rhs_Kn1 = dtn2t - dtn2;
 rhs_Kp1 = dtp2t - dtp2;
 rhs_G1  = dtH2t - dtH2;

 dtV2 = dtV2t;
 dtn2 = dtn2t;
 dtp2 = dtp2t;
 dtH2 = dtH2t;

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

 Ve = V + dt*( 1/6 * dtV1 + 2/3 * dtV2 + 1/6 * dtV3);
 ne = n + dt*( 1/6 * dtn1 + 2/3 * dtn2 + 1/6 * dtn3);
 pe = p + dt*( 1/6 * dtp1 + 2/3 * dtp2 + 1/6 * dtp3);
 He = H + dt*( 1/6 * dtH1 + 2/3 * dtH2 + 1/6 * dtH3);
 Ae = A + dt*( 1/6 *    H + 2/3 *   Hh + 1/6 *   He);

% [dtV3t,dtn3t,dtp3t,dtH3t,itNr3] = tdcalupdatec(Ve,ne,pe,Ae,He,dtVp,dtHp);
 [dtV3t,dtn3t,dtp3t,dtH3t] = tdcalupdates(Ve,ne,pe,Ae,He,dtVp,dtHp);

 rhs_F2  = dtV3t - dtV3;
 rhs_Kn2 = dtn3t - dtn3;
 rhs_Kp2 = dtp3t - dtp3;
 rhs_G2  = dtH3t - dtH3;

 dtV3 = dtV3t;
 dtn3 = dtn3t;
 dtp3 = dtp3t;
 dtH3 = dtH3t;

 rhs = [rhs_F1;rhs_Kn1;rhs_Kp1;rhs_G1;rhs_F2;rhs_Kn2;rhs_Kp2;rhs_G2];

 itNr = itNr+1;

 normUpdate = norm(rhs);

 display(['  Iter:',num2str(itNr),'; ','normRes:',num2str(normUpdate)]);

end

dtVt = (dtV1+4*dtV2+dtV3)/6;
dtnt = (dtn1+4*dtn2+dtn3)/6;
dtpt = (dtp1+4*dtp2+dtp3)/6;
dtHt = (dtH1+4*dtH2+dtH3)/6;
dtAt = (  H +4*Hh  +He  )/6;

Vn = V + dt* dtVt;
nn = n + dt* dtnt;
pn = p + dt* dtpt;
Hn = H + dt* dtHt;
An = A + dt* dtAt;


