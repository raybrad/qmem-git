%%%%%%%%%%%%%%%%%%%Parameter for generating meshgrids of whole box%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid_typ =1; %1 - rough ; 2 - fine(40); 3 - finer(60); 4 - finest(80)
metal_typ=0; %0 - no metal; 1 - plasmonic metal
noble_typ=2; %1 - Ag ; 2 - Au
if(grid_typ==1) 
dimensionL=40.0; %nm
latticeDL=2.0;  %nm
elseif(grid_typ==2) 
dimensionL=40.0; %nm
latticeDL=1.0;  %nm
elseif(grid_typ==3) 
dimensionL=40.0; %nm
latticeDL=0.6667d0;  %nm
elseif(grid_typ==4) 
dimensionL=40.0; %nm
latticeDL=0.5;  %nm
end

x_coor=[0:latticeDL:dimensionL];
y_coor=[0:latticeDL:dimensionL];
z_coor=[0:latticeDL:dimensionL];

plotx = false;
ploty = false;
plotz = false;

plotopt=2;  %1: plotting grid for the whole volume; 
            %2: plotting the first surface along each direction
            %3: plotting two end surfaces along each direction

gendata = true;  % true  = generate and output mesh data
                 % false = NOT do the work

meshfilename = 'meshGrids.mat';   %the output file will be in .MAT format

%%%%%%%%%%%%%%%%%%%Parameter for plasmonic metal and light%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (noble_typ==1) 
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

elseif (noble_typ==2) 
%Gold (From material database fitting with Drude+2 Lorentz model) 
PLANCKEV=4.13566733e-15;
epsilon_mt=3.559;
omega_p=8.812*2.0*pi/PLANCKEV;
gamma_p=0.0752*2.0*pi/PLANCKEV;

lepsr_1=2.912;
lomega_1=4.693*2.d0*pi/PLANCKEV;
lgamma_1=1.541*2.d0*pi/PLANCKEV;

lepsr_2=1.272;
lomega_2=3.112*2.d0*pi/PLANCKEV;
lgamma_2=0.525*2.d0*pi/PLANCKEV;
end

%light
dt  = 5.0e-17;
nsteps = 100;
epdf = 2.001347574e-15; %600nm epdf =wavelength/299.798
epdf2= 2.001347574e-15; %600nm
lightsource=1;
tzero=3.0e-15;
tlas=7.0e-16;
amplitude=6.0*1e6/(dimensionL*1e-9);
lightdirection='specialkzExBy';
Posz=35.0;
extfelc = 3; %zero voltage
%R
taskOpt =1;	%pure EM
%taskOpt=2;	%QM/EM
updateScheme=1; %implicit;
%updateScheme=2; %explicit;

nLead=0;
QMfeedback=0;

%%% three independent scaling parameters (lambda is changeable),default 1e-9
scl = struct('T',300,'Vt',2.5852,'lambda',1e-9,'ni',[],'s_D',[],'s_mu',[],...
    's_J',[],'tao',[],'s_E',[],'omega',[],'s_sigma',[],'s_Curr',[],...
    's_v',[],'s_A',[],'K',[]);

%%%%%%%%%%%%%%%%%%%Parameter for generating meshgrids of metal spherer%%%%%%%%%%%%%%%%%%%%%%%
Origin=[20.0,20.0,20.0];
Radius=10.0;
%%%%%%%parameter for output E Field%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outPutN=3;
cPosx=cell(outPutN,1);
cPosy=cell(outPutN,1);
cPosz=cell(outPutN+1,1);
for i=1:outPutN
    cPosx{i}=zeros(1,3);
    cPosy{i}=zeros(1,3);
    cPosz{i}=zeros(1,3);
   [cPosx{i}(1,1),cPosx{i}(1,2),cPosx{i}(1,3)]=RealToNodes([30+i,20,20],x_coor,y_coor,z_coor);
   [cPosy{i}(1,1),cPosy{i}(1,2),cPosy{i}(1,3)]=RealToNodes([20,30+i,20],x_coor,y_coor,z_coor);
   [cPosz{i}(1,1),cPosz{i}(1,2),cPosz{i}(1,3)]=RealToNodes([20,20,5*i],x_coor,y_coor,z_coor);
end
   [cPosz{outPutN+1}(1,1),cPosz{outPutN+1}(1,2),cPosz{outPutN+1}(1,3)]=RealToNodes([20,20,Posz],x_coor,y_coor,z_coor);
                                           

outputPosCom={cPosx{1},'Ex','11_xEx.dat';cPosx{2},'Ex','12_xEx.dat';cPosx{3},'Ex','13_xEx.dat';...
              cPosx{1},'By','11_xBy.dat';cPosx{2},'By','12_xBy.dat';cPosx{3},'By','13_xBy.dat';...
              cPosy{1},'Ex','11_yEx.dat';cPosy{2},'Ex','12_yEx.dat';cPosy{3},'Ex','13_yEx.dat';...
              cPosy{1},'By','11_yBy.dat';cPosy{2},'By','12_yBy.dat';cPosy{3},'By','13_yBy.dat';...
              cPosz{1},'Ex','5_zEx.dat' ;cPosz{2},'Ex','10_zEx.dat';cPosz{3},'Ex','15_zEx.dat';...
              cPosz{1},'By','5_zBy.dat' ;cPosz{2},'By','10_zBy.dat';cPosz{3},'By','15_zBy.dat';...
              cPosz{outPutN},'Ex','inc_Ex.dat';cPosz{outPutN},'By','inc_By.dat'};
tmpz=zRealTozNodes(20,z_coor);
outputPlane={'z',tmpz,'xyE2.dat'};
%%%%%%%%%%%%%%%%%%%Parameter for qm region%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[qmRegionX1,qmRegionY1,qmRegionZ1]=RealToNodes([30.0,18.0,18.0],x_coor,y_coor,z_coor);
[qmRegionX2,qmRegionY2,qmRegionZ2]=RealToNodes([34.0,22.0,22.0],x_coor,y_coor,z_coor);

[qmx1,qmy1,qmz1]=RealToNodes([30.0,18.0,18.0],x_coor,y_coor,z_coor);
[qmx2,qmy2,qmz2]=RealToNodes([34.0,22.0,22.0],x_coor,y_coor,z_coor);
%QM grid from lodestar
qkx =65 ; qky = 17; qkz = 17;
%R: set the boundary voltage to 0 now, for light
avoltage   = 0.0;
avoltageac = 0.0;
%%%%%%%
