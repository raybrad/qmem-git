clear global; clear;

global sigma0 sigma mu; %Conductivity and Permeability
sigma0 = 3.33e7^1;
mu = 4*pi*1e-7;
global epsilon epsilon_sd epsilon_mt epsilon_in epsilon_qm; %Dielectric constant
epsilon = 8.854e-12;
epsilon_sd = 11.9; %!undoped silicon substrate
epsilon_mt = 1;
epsilon_in = 3.9;  %!silicon dioxide
epsilon_qm = 1;
global ni; % intrinsic concentration
ni = 1e16;
global Ndp Nap; % carrier density of uniform doping
Ndp = 1.0e26;
Nap = -1.0e26;

global ehcrf;
ehcrf = 0.5;

global kb  %Boltzmann constant
kb = 1.3806488e-23;
global T  %Tempereture
T = 300;
global q  %Elementary charge
q = 1.6021892e-19;
global UT  %Boltzmann voltage
UT = kb*T/q;
global mun mup  %electrons and holes mobilities
mun=0.15;
mup=0.045;
global omega  %Frequency
global scl; 

global dt;
dt  = 1.0e-18;
global nsteps;
nsteps = 15;
global epdf;
epdf = 1.0e-15; 

global extfelc;
extfelc = 3;
global nedrelax;
nedrelax =2;
global irkod;
irkod  = 2;
global savefile;
savefile = 'dump/variables_';

global isehpgr;
isehpgr = 1;
global ehpgr;
ehpgr = 1.0e28;
global ehpkd ehpkf ehpkr;
%ehpkd = 5.0e6;
ehpkf = 5.0e12;
if mun < mup
  ehpkr = q*mun/(epsilon*epsilon_sd);
else
  ehpkr = q*mup/(epsilon*epsilon_sd); 
end

global ehpda;
ehpda = 1.6e-9;
global ebpkt;
ebpkt = exp(-q^2/(4*pi*epsilon*epsilon_sd*ehpda*kb*T));
ehpkd = 3*ehpkr*ebpkt/(4*pi*ehpda^3);

global ehdifc;
ehdifc = 1.0e-2;
%ehdifc = 0.0e-7;

global ehpkdns;

global fedstrgb;
fedstrgb = q^3/(8*pi*epsilon*epsilon_sd*kb^2*T^2);

%if isehpgr==1 
%  load('photocurrp');
%end



global kx ky kz;
global nodes links contacts;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
global Nnode Nlink Nsurf Nvolume;
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes;
global dirNodes eqnNodes;
global bndLinks  eqnLinks ;
global isBndNodes isDirNodes dcVolDirNodes acVolDirNodes; 
global doping;
global metalNodes eqnMetalNodes;
global semiNodes eqnSemiNodes dirSemiNodes isDirSemiNodes;
global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global isqmnodes;

global sQMlinks;
global currdlink;
global lqelinks rqelinks;
global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;
global ingvl ingvr;
global qxlinks;
global isqmvolm;


t_all = cputime;

%%%%% Load mesh data %%%%%%%%%%
loadPath = 'tdqmem.mat';  % path to load the mat file 
load(loadPath);

%%%%%%%%% geometry initialization %%%%%%%%
initGeometry;  

%%%%%%%%% scaling %%%%%%%%%%%%%%%
%scaling_const;
scaling_static; % scaling
scaling_dynamic(0);

% plotMesh;  % plot the mesh

%%%%%%%%% solver initialization (define material configuration) %%%%%%%%
initSolvertd;

qxlinks =[];

for kk=1:Nlink
    if isSQMlinks(kk) == 1 
      qxlinks  = [qxlinks;kk];
    end
end

%vbi = log(abs(Nap*Ndp/ni^2));
vbi = log(abs(Nap*Ndp/ni^2));
v3t = 1.0/scl.Vt;

  fp = fopen('ivcurve.dat','w');
  fclose(fp);

%%%%%%%%% static solution (can be skipped if no semiconductor material) %%%%%%%%%%%%%%%
%fav = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,1.0];
fav = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5];
Nf = length(fav);
  
for vav = 1:length(fav)

  savedump = ['dump',num2str(vav)];
  savefile = [savedump,'/variables_'];
  savetdcu = ['tdqmcurrd',num2str(vav),'.dat'];
  
  varcontv(fav(vav)/scl.Vt,0,vbi,0);

  display(['AppVolt:',num2str(fav(vav)),'; ','BuliV:',num2str(vbi*scl.Vt)]);

  initials = 'tdinis.mat';
  load(initials);

  ingvl=0;
  ingvr=0;

  [Vs,ns,ps,Xs] = staticSolutiond(V,n,p);
  V=Vs;
  n=ns;
  p=ps;
  X=Xs;
  save('tdinis.mat', 'V', 'n', 'p','X');

  [Ve,ne,pe,Xe,Ae,He] = statqmemx(dt,nsteps);
  load 'tdqmcurrd.dat';
  save(savetdcu,'tdqmcurrd','-ASCII');
  ! cd ./QMwork; ./remove.sh;


  currb=currentS([2,round(ky/2),round(kz/2)],Ve,ne,pe,1);

  fp = fopen('ivcurve.dat','a');
    fprintf(fp,'%f %f %f\n',fav(vav),currb*scl.s_Curr*1e18/4,currb*scl.s_Curr*1e9);
  fclose(fp);
  
end
t_all = cputime-t_all;
%figure; loglog(f,condNr);
