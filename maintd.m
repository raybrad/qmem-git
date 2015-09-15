%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Update History%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0 light source is added as a sheet of current density Js, which create plane wave field 
%0 Mur absorbing boundary is added to reduce reflection
%2 plasmonic metal is modelled using Two Lorentz poles and one drude pole.
%3 tdrelaxstep2.m use (row,col,val) technique to create sparse2 matrix
%3 in initGeometry sparse2 matrix creation is used for linkVolS nodeVolV
%4 tdbuilJacob.m prebuild,precondition and save  the Jacobian matrix,which is unchanged in time propagation. Semiconductor part is removed.
%4 tdsolverhs construct the right hand side matrix, and some simple vectorization technique is used. So tdrelaxstep2 is replaced.
%5 tdbuildRHSCoef fully vectorize the rhs coefficients,prebuild and save it. 
%5 tdsolerhs_vector use previous coefficient matrix to reconstruct the rhs matrix. So tdsolverhs is replaced. 
%6 in initGeometry the allocation of dirlinks is added to speed up
%6 in scaling_static change Vt from 2.5852e-2 to 2.5852,so to rescale a little bit to balance the num of Js by scl.J
%6 in initGeometry linkVolumes nodeVolumes are deleted and just used linkVolS and nodeVolV are ok,variables are changed in following sub
%7 in initGeometry nodeLinks construction is changed to sparse2 matrix creation to save time
%7 redundant global variables are removed(link(:,3) change also in initSolvertd
%7 unused subroutines is removed for clearance and there is no use to calculate main.m first,we could use V A H =0 as starting point.
%8 actually the linkVolumes nodeVolumes nodeLinks construction should be reserved as cells, since the use of find for sparse2 matrix in
%   later Jacob rhs construction is very time consuming. So now only reformulate the build up of linkvolumes and nodevolumes.
%9 fix a bug related to EsurfLinks and BsurfLinks,which is not allocated in specialkz case
%  plane wave generation is tested in no metal case,perfect planewave 

clear global; clear;
%plasmonic metal parameters
global sigma0  omega_p gamma_p; %Conductivity and Permeability
global mu;
global metal_typ;
global epsilon epsilon_mt epsilon_in epsilon_qm; %Dielectric constant
global lepsr_1 lomega_1 lgamma_1 lepsr_2 lomega_2 lgamma_2; %Lorentz pole


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
global nsteps;
global epdf epdf2;
global lightsource;
global tlas tzero;
global amplitude lightdirection Posz;  %k:the wave vector direction,E,B:polarization  
global outputPosCom outputPlane;

global extfelc;
global nedrelax;
global nLead;
global QMfeedback;
global irkod;
irkod  = 2;
global savefile;
savefile = 'dump/variables_';

global x_coor y_coor z_coor dimensionL latticeDL;
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
global qmRegionX1 qmRegionX2 qmRegionY1 qmRegionY2 qmRegionZ1 qmRegionZ2;
global ingvl ingvr ingvt;
global qxlinks;
global isqmvolm;
global qkx qky qkz;
global avoltage avoltageac;

t_all = cputime;
%%%load Parameters%%%%%%%%%%%%%
inputParaDefault;
inputParaUser;

%loadPath='inputPara.m';
%load(loadPath);
%%%generate Mesh Grids%%%%%%%%%
genMeshGrids(x_coor,y_coor,z_coor,plotx,ploty,plotz,plotopt,gendata,meshfilename);

%%%%% Load mesh data %%%%%%%%%%
load(meshfilename);

%%%%%%%%% geometry initialization %%%%%%%%
initGeometry;  

%%%%%%%%% scaling %%%%%%%%%%%%%%%
%scaling_const;
scaling_static; % scaling
scaling_dynamic(0);


%%%%%generate mesh for metal sphere%%%%%%%%%%%%%%%%%%
if (metal_typ== 1) 
metalNodes=defSphereNodes(Origin,scl.lambda,Radius,nodes);
elseif (metal_typ==0) 
metalNodes=[];
else
    fprintf('other Metal_typ not defined now \n');
end
%savefilename='metalNodes.mat';
%save(savefilename, 'metalNodes','Nnode','kx','ky','kz','x_coor','y_coor','z_coor');
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
