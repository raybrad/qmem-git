function [V,A,H] = tdsimulSolution(dt,nsteps)

display(['Start TD dynamic solution']);

global nedrelax;
global irkod;
global savefile;
global Nlink;
global Nnode;
global extfelc;
global dirNodes dcVolDirNodes acVolDirNodes;
global epdf epdf2;

global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global isqmnodes;

global JLinks J_amp;
global lightsource tlas tzero;
global outputPosCom outputPlane;
%global Jacob colind Lmatrix Umatrix;
%global JacobLU;
%global Fc11 Fc12 Fc21 Fc22 Fc30 Fn1matrix Fn2matrix Flkmatrix;
%global Gc11 Gajlkmatrix11 Gc12 Gc20 Gn1matrix Gc31 Gajlk_n1matrix Gc32;
%global Gn2matrix Gc41 Gajlk_n2matrix Gc42 Gc51 Gc52 Gc53 Gc54;
%%%%%%% initial guess  %%%%%%%%%%%%%
%initials = 'tdinis.mat';

%load(initials);
V=zeros(Nnode,1);
H=zeros(Nlink,1);
A=zeros(Nlink,1);

interval=10;

mJ_0=zeros(Nlink,1);
mJ_1=zeros(Nlink,1);
mJ_2=zeros(Nlink,1);
mJ_1p=zeros(Nlink,1);
mJ_2p=zeros(Nlink,1);
Efield_p=zeros(Nlink,1);

%R
Js=zeros(Nlink,1);
dtVp = zeros(Nnode,1);
dtHp = zeros(Nlink,1);

timec=0;

switch extfelc
   case 1
     dtVp(dirNodes) = acVolDirNodes*exp(- timec/epdf)/epdf;
   case 2
     dtVp(dirNodes) = pi*acVolDirNodes*sin(2*pi*timec/epdf)/epdf; 
   case 3
     dtVp(dirNodes) = 0.0;
   case 4
     dtVp(dirNodes) = -acVolDirNodes*exp((timec-nsteps*dt)/epdf)/epdf;
   case 5
     dtVp(dirNodes) = -2*(timec-0.5*nsteps*dt)*acVolDirNodes*exp(-((timec-0.5*nsteps*dt)/epdf)^2)/epdf^2;
   otherwise
     error('undefined bias profile');
end

% a current density sheet is added to simulate the light source
switch lightsource
 case 1
 Js(JLinks)=J_amp*sin(2*pi*timec/epdf2);
 case 2
 Js(JLinks)=J_amp*sin(2*pi*timec/epdf2)*exp(-((timec-tzero)/tlas)^2.0);   %omega= 2*pi/1 ,E = hbar *omega =0.6582*2*3.14/1=4.1356eV ~ 300nm
 case 3
 if (timec<=7*tlas) 
 Js(JLinks)=J_amp*exp(-((timec-tzero)/tlas)^2.0);   %omega= 2*pi/1 ,E = hbar *omega =0.6582*2*3.14/1=4.1356eV ~ 300nm
 else
 Js(JLinks)=0.0;
 end
 otherwise
 error('undefined light source');
end

%tdcalupdatec.m updates dtV,dtH,dtn,dtp 
%(more specificly,Newton's method to get dtV dtH.. dtn dtp are calculated directly usiing continuity equation)
%tdrelaxstep2.m and others update  V n p A H (dV,dn ,dp,dA,dH are in Newton)(dtV,dtH,dtn,dtp are calculated by  numerical differentiation)
if (nedrelax==1)
[dtVp,dtHp] = tdcalupdatec(V,A,H,Js,dtVp,dtHp);
elseif (nedrelax==2)
[dtVp,dtHp] = tdcalupdatecc(V,A,H,Js,dtVp,dtHp);
end

%savefilename = [savefile,num2str(0),'.mat'];
%save(savefilename, 'V', 'A', 'H','mJ_0','mJ_1','mJ_2','mJ_1p','mJ_2p','dtVp','dtHp','Js');
% save(savefilename, 'V','A', 'H','dtVp','dtHp','Js');

if nedrelax==2
  fp = fopen('tdqmcurrd.dat','w');  % clear content of dat & create file
  fclose(fp);
end

ntp  = 1; 
tic;
%if  (exist('JacobMatrix.mat') == 0 )
%tdbuildJacob(dt);
%else
%display('JacobMatrix.mat already exist');
%end
%load 'JacobMatrix.mat';
if (nedrelax == 1)
%if  (exist('JacobLU.mat') == 0 )
tdbuildJacob(dt);
%else
%display('JacobLU.mat already exist');
%end

%load 'JacobLU.mat';

%if  (exist('rhsFCoefMatrix.mat') == 0 || exist('rhsGCoefMatrix.mat') == 0 )
tdbuildRHSCoef(dt);
%else
%display('rhsF/GCoefMatrix.mat already exist');
%end
elseif (nedrelax==2)
tdbuildJacobc(dt);
tdbuildRHSCoefc(dt);
end
%load 'rhsFCoefMatrix.mat';
%load 'rhsGCoefMatrix.mat';

while ntp < nsteps + 1
   
 display(['Time step:',num2str(ntp)]);

 if nedrelax==1
   [Vu,Au,Hu,mJ_0u,mJ_1u,mJ_2u,mJ_1uu,mJ_2uu,Efield_uu,dtVu,dtHu,Js]=tdupdatenm(V,A,H,mJ_0,mJ_1,mJ_2,mJ_1p,mJ_2p,Efield_p,dtVp,dtHp,dt,ntp);
 elseif nedrelax==2
   [Vu,Au,Hu,mJ_0u,mJ_1u,mJ_2u,mJ_1uu,mJ_2uu,Efield_uu,dtVu,dtHu,Js]=tdupdatenmc(V,A,H,mJ_0,mJ_1,mJ_2,mJ_1p,mJ_2p,Efield_p,dtVp,dtHp,dt,ntp);
 end

 dtVp  = dtVu;
 dtHp  = dtHu;
 
 V = Vu;
 A = Au;
 H = Hu;

 mJ_0= mJ_0u;	%J_n+1 = J_n
 mJ_1= mJ_1u;
 mJ_2= mJ_2u;
 mJ_1p = mJ_1uu; %J_n = J_n-1
 mJ_2p = mJ_2uu;
 Efield_p=Efield_uu;

%if (ntp>10 &&  mod(ntp,interval) == 0 )
% savefilename = [savefile,num2str(ntp),'.mat'];
% save(savefilename, 'V','A', 'H','dtVp','dtHp','Js');
%elseif (ntp<=10)
% savefilename = [savefile,num2str(ntp),'.mat'];
% save(savefilename, 'V','A', 'H','Js');
%end
outputField(ntp,outputPosCom,outputPlane,V,A,H);

ntp = ntp + 1;
    
end    
toc;   

display(['End TD dynamic solution']);
