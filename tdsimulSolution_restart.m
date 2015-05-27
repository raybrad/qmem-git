function [V,n,p,A,H] = tdsimulSolution_restart(dt,nsteps)

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
%%%%%%% initial guess  %%%%%%%%%%%%%
%V = Vs;
%n = ns;
%p = ps;
%A = zeros(Nlink,1);
%H = zeros(Nlink,1);

%[A,H] = tdstatica(V,n,p);

%initials = 'tdinis.mat';
initials='dump/variables_550.mat';

load(initials);

mJ=zeros(Nlink,1);

Vpp = V;
npp = n;
ppp = p;
App = A;
Hpp = H;
%R
Js=zeros(Nlink,1);
dtVp = zeros(Nnode,1);
dtnp = zeros(Nnode,1);
dtpp = zeros(Nnode,1);
dtHp = zeros(Nlink,1);

dtVp2= zeros(Nnode,1);
dtnp2= zeros(Nnode,1);
dtpp2= zeros(Nnode,1);
dtHp2= zeros(Nlink,1);

%timec=0;

%switch extfelc
%   case 1
%     dtVp(dirNodes) = acVolDirNodes*exp(- timec/epdf)/epdf;
%   case 2
%     dtVp(dirNodes) = pi*acVolDirNodes*sin(2*pi*timec/epdf)/epdf; 
%   case 3
%     dtVp(dirNodes) = 0.0;
%   case 4
%     dtVp(dirNodes) = -acVolDirNodes*exp((timec-nsteps*dt)/epdf)/epdf;
%   case 5
%     dtVp(dirNodes) = -2*(timec-0.5*nsteps*dt)*acVolDirNodes*exp(-((timec-0.5*nsteps*dt)/epdf)^2)/epdf^2;
%   otherwise
%     error('undefined bias profile');
%end

% a current density sheet is added to simulate the light source
 Js(JLinks)=J_amp*sin(2*pi*timec/epdf2);   %omega= 2*pi/1 ,E = hbar *omega =0.6582*2*3.14/1=4.1356eV ~ 300nm

%tdcalupdatec.m updates dtV,dtH,dtn,dtp 
%(more specificly,Newton's method to get dtV dtH.. dtn dtp are calculated directly usiing continuity equation)
%tdrelaxstep2.m and others update  V n p A H (dV,dn ,dp,dA,dH are in Newton)(dtV,dtH,dtn,dtp are calculated by  numerical differentiation)
%[dtVp,dtnp,dtpp,dtHp,itNr1] = tdcalupdatec(V,n,p,A,H,mJ,Js,dtVp,dtHp);

% dtVp2 = dtVp;
% dtnp2 = dtnp;
% dtpp2 = dtpp;
% dtHp2 = dtHp;

%savefilename = [savefile,num2str(0),'.mat'];
%save(savefilename, 'V', 'n', 'p', 'A', 'H','mJ','dtVp','dtnp','dtpp','dtHp','Js');

if nedrelax==2
  fp = fopen('tdqmcurrd.dat','w');  % clear content of dat & create file
  fclose(fp);
end

ntp  = 551; 
tic;

while ntp < nsteps + 1
   
 display(['Time step:',num2str(ntp)]);

 if nedrelax==0
    [Vu,nu,pu,Au,Hu,dtVu,dtnu,dtpu,dtHu] = tdupdaterk(V,n,p,A,H,dtVp,dtnp,dtpp,dtHp,dt,ntp);
 elseif nedrelax==1
    [Vu,nu,pu,Au,Hu,mJu,dtVu,dtnu,dtpu,dtHu,Js] = tdupdatenm(V,n,p,A,H,mJ,Vpp,npp,ppp,App,Hpp,dtVp,dtnp,dtpp,dtHp,dtVp2,dtnp2,dtpp2,dtHp2,dt,ntp);
 elseif nedrelax==2
    [Vu,nu,pu,Au,Hu,mJu,dtVu,dtnu,dtpu,dtHu,Js] = tdupdatenmc(V,n,p,A,H,mJ,Vpp,npp,ppp,App,Hpp,dtVp,dtnp,dtpp,dtHp,dtVp2,dtnp2,dtpp2,dtHp2,dt,ntp);
 end

 dtVp2 = dtVp;
 dtnp2 = dtnp;
 dtpp2 = dtpp;
 dtHp2 = dtHp;

 dtVp  = dtVu;
 dtnp  = dtnu;
 dtpp  = dtpu;
 dtHp  = dtHu;

 Vpp = V;
 npp = n;
 ppp = p;
 App = A;
 Hpp = H;
 
 V = Vu;
 n = nu;
 p = pu;
 A = Au;
 H = Hu;
 mJ= mJu;
    
 savefilename = [savefile,num2str(ntp),'.mat'];
 save(savefilename, 'V', 'n', 'p', 'A', 'H','mJ','dtVp','dtnp','dtpp','dtHp','Js');

    ntp = ntp + 1;
    
end    
toc;   

display(['End TD dynamic solution']);
