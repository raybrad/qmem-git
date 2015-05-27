function [Vn,nn,pn] = statupdatenmc(V,n,p,A,H,dt,ntp)

global epdf;
global extfelc;
global dirNodes dcVolDirNodes acVolDirNodes;
global nsteps;
global irkod;

global currdlink;
global bqelinks tqelinks;
global lqelinks rqelinks;
global qxlinks;
global scl;

timec = (ntp-1)*dt;
timeh = timec + 0.5*dt;
timee = timec + dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pass Vs to LODESTAR
%
% call lodestar
qmeffcd = scl.ni*effchgden(H);
qmv  = V*scl.Vt;
qmdt = dt*scl.tao*1e15;
currtime = timec*scl.tao*1e15;
[currdl,currdr,currdt]=QM_Calc(qmv,qmeffcd,currtime,qmdt);
currdlink(qxlinks,4)=(currdl+currdr)/(2*scl.s_J);
%z direction?

fp = fopen('tdqmcurrd.dat','a');
fprintf(fp,'%f %.15f %.15f %.15f %.15f\n',currtime,currdl,currdr,currdt,(currdl+currdr)/2);
fclose(fp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Vn,nn,pn] = staticSolutionc(V,n,p,currdlink);
