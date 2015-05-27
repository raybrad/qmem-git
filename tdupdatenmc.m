function [Vn,An,Hn,mJ_0n,mJ_1n,mJ_2n,mJ_1nn,mJ_2nn,Efield_nn,dtVu,dtHu,Js] = tdupdatenmc(V,A,H,mJ_0,mJ_1,mJ_2,mJ_1p,mJ_2p,Efield_p,dtVp,dtHp,dt,ntp)

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

%%@
switch lightsource
 case 1
 Js(JLinks)=J_amp*sin(2*pi*timee/epdf2);
 case 2
 Js(JLinks)=J_amp*sin(2*pi*timee/epdf2)*exp(-((timee-tzero)/tlas)^2.0);   %omega= 2*pi/1 ,E = hbar *omega =0.6582*2*3.14/1=4.1356eV ~ 300nm
 case 3
 if (timee<=7*tlas) 
 Js(JLinks)=J_amp*exp(-((timee-tzero)/tlas)^2.0);   %omega= 2*pi/1 ,E = hbar *omega =0.6582*2*3.14/1=4.1356eV ~ 300nm
 else
 Js(JLinks)=0.0;
 end
 otherwise
 error('undefined light source');
end
%%%
[Vn,An,Hn,mJ_0n,mJ_1n,mJ_2n,mJ_1nn,mJ_2nn,Efield_nn,dtVu,dtHu] =  tdsolverhs_vectorc(V,A,H,Js,mJ_0,mJ_1,mJ_2,mJ_1p,mJ_2p,Efield_p,dtV,dtH,dt);

%[Vn,An,Hn,mJ_0n,mJ_1n,mJ_2n,mJ_1nn,mJ_2nn,Efield_nn,dtVu,dtHu] =  tdrelaxpsteps2c(V,A,H,Js,mJ_0,mJ_1,mJ_2,mJ_1p,mJ_2p,Efield_p,dtV,dtH,dt);
%%%[dtVu,dtnu,dtpu,dtHu,itNr] = tdcalupdatecc(Vn,nn,pn,An,Hn,dtVp,dtHp);
