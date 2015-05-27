clear global; clear;

%plasmonic metal parameters
global sigma0  omega_p gamma_p; %Conductivity and Permeability Drude pole
global mu;
global epsilon epsilon_mt epsilon_in epsilon_qm; %Dielectric constant
global lepsr_1 lomega_1 lgamma_1 lepsr_2 lomega_2 lgamma_2; %Lorentz pole

%Silver (John and Christy:400nm-800nm Drude model)
%epsilon_mt= 4.07666;
%omega_p   = 1.40056E16; 
%gamma_p   = 4.21755E13;

%Silver (From material database fitting with Drude+2 Lorentz model) 
%suitable over 1.24 - 5 eV

PLANCKEV=4.13566733e-15;
epsilon_mt=3.189;
omega_p=9.183*2.0*pi/PLANCKEV;
gamma_p=0.0179*2.0*pi/PLANCKEV;

lepsr_1=0.4323;
lomega_1=4.668*2.0*pi/PLANCKEV;
lgamma_1=0.207*2.0*pi/PLANCKEV;

lepsr_1=0.2237;
lomega_1=4.240*2.0*pi/PLANCKEV;
lgamma_1=0.186*2.0*pi/PLANCKEV;

mu = 4*pi*1e-7;
epsilon = 8.854e-12; %vaccum(absolute)

epsilon_in = 1.0;    %insulator
epsilon_qm = 1;      %qm region 

%not used,just for  comparison 
sigma0=epsilon*omega_p*omega_p/gamma_p;

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
global omega  %Frequency
global scl; 

global dt;
dt  = 5.0e-17; %0.01fs
global nsteps;
nsteps = 10000;
global epdf epdf2;
epdf = 2.001347574e-15; %600nm epdf =wavelength/299.798
epdf2= 2.001347574e-15; %600nm
global lightsource;
lightsource=2;
global tlas tzero;
tzero=10.0e-15;
tlas=7.0e-16;
global lightdirection Posz;  %k:the wave vector direction,E,B:polarization  
lightdirection='specialkzExBy';
Posz=35.0;

global extfelc;
extfelc = 3; %zero voltage
global nedrelax;
%R
nedrelax =1;	%pure EM
%nedrelax=2;	%QM/EM
global irkod;
irkod  = 2;
global savefile;
savefile = 'dump/variables_';

global kx ky kz;
global nodes links contacts;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolS nodeVolV;
global Nnode Nlink Nsurf Nvolume;
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes;
global dirNodes eqnNodes;
global bndLinks  eqnLinks ;
global isBndNodes isDirNodes dcVolDirNodes acVolDirNodes; %#ok<*NUSED>
global doping;
global metalNodes eqnMetalNodes;
global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global isqmnodes;

global sQMlinks;
global currdlink;
global lqelinks rqelinks bqelinks tqelinks;
global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;
global ingvl ingvr ingvt;
global qxlinks;
global isqmvolm;


t_all = cputime;

%%%%% Load mesh data %%%%%%%%%%
loadPath = 'tdmetal.mat';  % path to load the mat file 
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
tqelinks=[];
qxlinks =[];

for kk=1:Nlink
    if isSQMlinks(kk) == 1 
      lqelinks = [lqelinks;kk];
    elseif isSQMlinks(kk) == 2
      rqelinks = [rqelinks;kk];
    elseif isSQMlinks(kk) == 6
      tqelinks = [tqelinks;kk];
    elseif isSQMlinks(kk) == 7
      qxlinks  = [qxlinks;kk];
    end
end

%%%%%%%%% static solution (can be skipped if no semiconductor material) %%%%%%%%%%%%%%%

ingvl=0;
ingvr=0;
ingvt=0;

%%stop;


%%%%%%%%% static solution (can be skipped if no semiconductor material) %%%%%%%%%%%%%%%
[Vs] = staticSolution;  %to get V n p


[A,H] = tdstatica(Vs);  %to get A H~par A/par t

V=Vs;

save('tdinis.mat', 'V', 'A', 'H');

t_all = cputime-t_all;
%figure; loglog(f,condNr);
