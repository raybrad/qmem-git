function [Vn,nn,pn,An,Hn,mJn,dtVu,dtnu,dtpu,dtHu,Js] = tdupdatenm(V,n,p,A,H,mJ,Vpp,npp,ppp,App,Hpp,dtVp,dtnp,dtpp,dtHp,dtVp2,dtnp2,dtpp2,dtHp2,dt,ntp)

global epdf epdf2;
global extfelc;
global dirNodes dcVolDirNodes acVolDirNodes;
global nsteps;
global irkod;
global lightsource tlas tzero;

global currdlink;
global bqelinks tqelinks;
global lqelinks rqelinks;
global scl;
global qxlinks;
%%
global Nlink JLinks J_amp;
%%
Js=zeros(Nlink,1);
timec = (ntp-1)*dt; %
timeh = timec + 0.5*dt;
timee = timec + dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pass Vs to LODESTAR
%
% call lodestar
qmeffcd = scl.ni*effchgden(H);     %integral_S div dA/dt -> effective charge density
qmv  = V*scl.Vt;
qmdt = dt*scl.tao*1e15;
currtime = timec*scl.tao*1e15;
[currdl,currdr,currdt]=QM_Calc(qmv,qmeffcd,currtime,qmdt);
currdlink(lqelinks,4)=currdl/scl.s_J;
currdlink(rqelinks,4)=currdr/scl.s_J;
currdlink(tqelinks,4)=currdt/scl.s_J;
currdlink(qxlinks,4)=(currdl+currdr)/(2*scl.s_J);
%qzlinks


fp = fopen('tdqmcurrd.dat','a');
fprintf(fp,'%f %.15f %.15f  %.15f %.15f\n',currtime,currdl,currdr,currdt,(currdl+currdr)/2);
fclose(fp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 dtV = dtVp;
 dtn = dtnp;
 dtp = dtpp;
 dtH = dtHp;
switch extfelc
   case 1
     dtVp(dirNodes) = acVolDirNodes*exp(- timeh/epdf)/epdf;
   case 2
     dtVp(dirNodes) = pi*acVolDirNodes*sin(2*pi*timeh/epdf)/epdf;
   case 3
     dtVp(dirNodes) = 0.0;
   case 4
     dtVp(dirNodes) = -acVolDirNodes*exp((timeh-nsteps*dt)/epdf)/epdf;
   case 5
     dtVp(dirNodes) = -2*(timeh-0.5*nsteps*dt)*acVolDirNodes*exp(-((timeh-0.5*nsteps*dt)/epdf)^2)/epdf^2;
   otherwise
     error('undefined bias profile');
end

dtV2 = dtVp;

switch extfelc
   case 1
     dtVp(dirNodes) = acVolDirNodes*exp(- timee/epdf)/epdf;
   case 2
     dtVp(dirNodes) = pi*acVolDirNodes*sin(2*pi*timee/epdf)/epdf;
   case 3
     dtVp(dirNodes) = 0.0;
   case 4
     dtVp(dirNodes) = -acVolDirNodes*exp((timee-nsteps*dt)/epdf)/epdf;
   case 5
     dtVp(dirNodes) = -2*(timee-0.5*nsteps*dt)*acVolDirNodes*exp(-((timee-0.5*nsteps*dt)/epdf)^2)/epdf^2;
   otherwise
     error('undefined bias profile');
end

dtV3 = dtVp;
dtVt = (dtV+4*dtV2+dtV3)/6;
Vn = V + dt * dtVt;
dtV(dirNodes)  = dtVp(dirNodes);
V= Vn - dt * dtV;

%%%@
switch lightsource
 case 1
 Js(JLinks)=J_amp*sin(2*pi*timee/epdf2);
 case 2
 Js(JLinks)=J_amp*sin(2*pi*timee/epdf2)*exp(-((timee-tzero)/tlas)^2.0);   %omega= 2*pi/1 ,E = hbar *omega =0.6582*2*3.14/1=4.1356eV ~ 300nm
 otherwise
 error('undefined light source');
end
%%%
switch irkod
   case 1  
     [Vn,nn,pn,An,Hn] = tdrelaxpstepsc(V,n,p,A,H,dtV,dtn,dtp,dtH,dt);
   case 2  %
     [Vn,nn,pn,An,Hn,mJn,dtVu,dtnu,dtpu,dtHu] = tdrelaxpsteps2c(V,n,p,A,H,Js,mJ,dtV,dtn,dtp,dtH,dt);
   case 4
     [Vn,nn,pn,An,Hn] = tdrelaxpsteps4c(V,n,p,A,H,Vpp,npp,ppp,App,Hpp,dtV,dtn,dtp,dtH,dtVp2,dtnp2,dtpp2,dtHp2,dt);
     dtVp  = 3 * (Vn - V)/dt - 4*dtV - dtVp2;
     dtHp  = 3 * (Hn - H)/dt - 4*dtH - dtHp2;
     dtVp(dirNodes)  = dtV(dirNodes);
   otherwise
     error('undefined irkod');
end

%[dtVu,dtnu,dtpu,dtHu,itNr] = tdcalupdatecc(Vn,nn,pn,An,Hn,dtVp,dtHp);
