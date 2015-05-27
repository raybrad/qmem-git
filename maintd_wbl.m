% close all;
clear global; clear;

global sigma0 sigma mu; %Conductivity and Permeability
sigma0 = 3.33e7^1;
mu = 4*pi*1e-7;
global epsilon epsilon_sd epsilon_mt epsilon_in epsilon_qm; %Dielectric constant
epsilon = 8.854e-12;
epsilon_sd = 11.9;
epsilon_mt = 1;
epsilon_in = 3.9;
epsilon_qm = 1;
global ni; % intrinsic concentration
ni = 1e16;
global Ndp Nap; % carrier density of uniform doping
Ndp = 1e24;
Nap = 1e24;

global ehcrf;
ehcrf = 0.5;

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
dt  = 1.5e-17;
global nsteps;
nsteps = 8000;
global epdf;
epdf = 5.0e-15; 

global extfelc;
extfelc = 2;
global nedrelax;
nedrelax =2;
global irkod;
irkod  = 2;
global savefile;
savefile = 'dump/variables_';

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
global lqelinks rqelinks;
global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;
global ingvl ingvr;
global qxlinks;
global isqmvolm;


t_all = cputime;

%%%%% Load mesh data %%%%%%%%%%
loadPath = 'tdqmem.mat';  % path to load the mat file 
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

lqelinks=[];
rqelinks=[];
qxlinks =[];

for kk=1:Nlink
    if isSQMlinks(kk) == 1 
      lqelinks = [lqelinks;kk];
    elseif isSQMlinks(kk) == 2
      rqelinks = [rqelinks;kk];
    elseif isSQMlinks(kk) == 7
      qxlinks  = [qxlinks;kk];
    end
end

%%%%%%%%% static solution (can be skipped if no semiconductor material) %%%%%%%%%%%%%%%

ingvl=0;
ingvr=0;

[Ve,ne,pe,Ae,He] = tdsimulSolution(dt,nsteps);

[VT,nT,pT,AT,HT] = tdvara(nsteps);

t_all = cputime-t_all;
%figure; loglog(f,condNr);
