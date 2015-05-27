function genMetalNodes
global kx ky kz;
load 'tdmetal.mat';
nodes = nodes(:,2:4);
kx = length(unique(nodes(:,1)))-1;
ky = length(unique(nodes(:,2)))-1;
kz = length(unique(nodes(:,3)))-1;
lambda=1e-9;


x_coor=[0:1:9  9.5:0.5:30.5 31:1:40];
y_coor=[0:1:9  9.5:0.5:30.5 31:1:40];
z_coor=[0:1:9  9.5:0.5:30.5 31:1:40];

Radius=10.0;
metalNodes=defSphereNodes([20.0,20.0,20.0],lambda,Radius,nodes);

savefilename='metalNodes2.mat';
save(savefilename, 'metalNodes');
%plotStructureIso;
