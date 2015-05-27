% close all;
clear global; clear;

global sigma0 sigma mu; %Conductivity and Permeability
sigma0 = 3.33e7^1;
mu = 4*pi*1e-7;
global epsilon epsilon_sd epsilon_mt epsilon_in; %Dielectric constant
epsilon = 8.854e-12;
epsilon_sd = 11.9;
epsilon_mt = 1;
epsilon_in = 3.9;
global ni; % intrinsic concentration
ni = 1e16;
global Ndp Nap; % carrier density of uniform doping
Ndp = 1e24;
Nap = 1e24;

global ehcrf;
ehcrf = 1.0;

global kb  %Boltzmann constant
kb = 1.38e-23;
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
dt  = 2.0e-14;
global nsteps;
nsteps = 50;
global epdf;
epdf = 1.0e-12; 

global extfelc;
extfelc = 3;
global nedrelax;
nedrelax =1;
global irkod;
irkod  = 1;
global savefile;
savefile = './compt2/src/dump/variables_';

global kx ky kz;
global nodes links contacts;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
global Nnode Nlink Nsurf Nvolume;
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes;
global dirNodes eqnNodes;
global bndLinks  eqnLinks ;
global isBndNodes isDirNodes dcVolDirNodes acVolDirNodes; %#ok<*NUSED>
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
global isqmvolm;
global epsilon_qm;
epsilon_qm = 1;


t_all = cputime;

%%%%% Load mesh data %%%%%%%%%%
loadPath = './cntal_5kblock/tdqmem.mat';  % path to load the mat file 
%loadPath = 'meshq4.mat';
% loadPath = 'mesh_metalplug.mat';
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

nsteps=10000;

[VT,nT,pT,AT,HT]=tdvara(nsteps);
dembv=tdembv(VT);
[currentT1,currddT1]=tdcurrentf([3,6,6],1,nsteps);
[currentT2,currddT2]=tdcurrentfqmd([8,6,6],1,nsteps);
%[currentT3,currddT3]=tdcurrentf([3,6,6],1,nsteps);
%[currentT4,currddT4]=tdcurrentf([kx-2,6,6],1,nsteps);
%load 'QMwork/curr.data';
load './compt2/src/tdqmcurrd.dat';
%ansv1=[currentT1,currentT3,currentT2];
%ansv2=[VT(901,:)',currentT1,currentT2,currentT3];
ansv3=[scl.Vt*VT(901,:)',scl.Vt*dembv,currentT1*scl.s_Curr*1e9,currentT2*scl.s_Curr*1e9,(currentT1-currentT2)*scl.s_Curr*1e9,tdqmcurrd(1:nsteps,2)*0.75*0.75*1e-9];

save('currqmem4.dat','ansv3','-ASCII');

%[Ve,ne,pe,Ae,He] = tdsimulSolution(dt,nsteps);

%[VT,nT,pT,AT,HT] = tdvara(nsteps);

t_all = cputime-t_all;
%figure; loglog(f,condNr);
