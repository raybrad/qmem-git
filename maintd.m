%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Update History%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0 light source is added as a sheet of current density Js, which create plane wave field 
%0 Mur absorbing boundary is added to reduce reflection
%2 plasmonic metal is modelled using Two Lorentz poles and one drude pole.
%3 tdrelaxstep2.m use (row,col,val) technique to create sparse matrix
%3 in initGeometry sparse matrix creation is used for linkVolS nodeVolV
%4 tdbuilJacob.m prebuild,precondition and save  the Jacobian matrix,which is unchanged in time propagation. Semiconductor part is removed.
%4 tdsolverhs construct the right hand side matrix, and some simple vectorization technique is used. So tdrelaxstep2 is replaced.
%5 tdbuildRHSCoef fully vectorize the rhs coefficients,prebuild and save it. 
%5 tdsolerhs_vector use previous coefficient matrix to reconstruct the rhs matrix. So tdsolverhs is replaced. 
%6 in initGeometry the allocation of dirlinks is added to speed up
%6 in scaling_static change Vt from 2.5852e-2 to 2.5852,so to rescale a little bit to balance the num of Js by scl.J
%6 in initGeometry linkVolumes nodeVolumes are deleted and just used linkVolS and nodeVolV are ok,variables are changed in following sub
%7 in initGeometry nodeLinks construction is changed to sparse matrix creation to save time
%7 redundant global variables are removed(link(:,3) change also in initSolvertd
%7 unused subroutines is removed for clearance and there is no use to calculate main.m first,we could use V A H =0 as starting point.
%8 actually the linkVolumes nodeVolumes nodeLinks construction should be reserved as cells, since the use of find for sparse matrix in
%   later Jacob rhs construction is very time consuming. So now only reformulate the build up of linkvolumes and nodevolumes.
%9 use tim davis's sparse2 method to replace sparse construction. and multiple linear solvers are tested,but lu remains best of them
%10 reformulate the tdcalupdatec and tdbuildJaocb by part vectorlization,but still not fast,now is a bottleneck for finer grids

clear global; clear;

%plasmonic metal parameters
global sigma0  omega_p gamma_p; %Conductivity and Permeability
global mu;
global epsilon epsilon_mt epsilon_in epsilon_qm; %Dielectric constant
global lepsr_1 lomega_1 lgamma_1 lepsr_2 lomega_2 lgamma_2; %Lorentz pole
%Silver (From material database fitting with Drude+2 Lorentz model) 
%suitable over 1.24 - 5 eV
PLANCKEV=4.13566733e-15;
epsilon_mt=3.189;
omega_p=9.183*2.0*pi/PLANCKEV;
gamma_p=0.0179*2.0*pi/PLANCKEV;

lepsr_1=0.4323;
lomega_1=4.668*2.d0*pi/PLANCKEV;
lgamma_1=0.207*2.d0*pi/PLANCKEV;

lepsr_2=0.2237;
lomega_2=4.240*2.d0*pi/PLANCKEV;
lgamma_2=0.186*2.d0*pi/PLANCKEV;

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
dt  = 5.0e-17;
global nsteps;
nsteps = 1000;
global epdf epdf2;
epdf = 2.001347574e-15; %600nm epdf =wavelength/299.798
epdf2= 2.001347574e-15; %600nm
global lightsource;
lightsource=2;
global tlas tzero;
tzero=3.0e-15;
tlas=7.0e-16;
global lightdirection Posz;  %k:the wave vector direction,E,B:polarization  
lightdirection='specialkzExBy';
Posz=35;

global extfelc;
extfelc = 3; %zero voltage
global nedrelax;
%R
nedrelax =1;	%pure EM
%nedrelax=2;	%QM/EM
global nLead;
nLead=0;
global irkod;
irkod  = 2;
global savefile;
savefile = 'dump/variables_';

global kx ky kz;
global nodes links contacts;
global nodeLinks linkSurfs  surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
global Nnode Nlink Nsurf Nvolume;
global nodeV linkL linkS dlinkL nodeM linkM volumeM;
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

[Ve,Ae,He] = tdsimulSolution(dt,nsteps);

t_all = cputime-t_all;
