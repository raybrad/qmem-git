
%clc;
%
%x_coor = [0.0:1.0:15.0];
%y_coor = [0.0:1.0:15.0];
%z_coor = [0.0:1.0:15.0];
% 81
x_coor=[0:1.0:40];
y_coor=[0:1.0:40];
z_coor=[0:1.0:40];

plotx = false;
ploty = false;
plotz = false;

plotopt=2;  %1: plotting grid for the whole volume; 
            %2: plotting the first surface along each direction
            %3: plotting two end surfaces along each direction

gendata = true;  % true  = generate and output mesh data
                 % false = NOT do the work

savefilename = 'tdmetal.mat';   %the output file will be in .MAT format

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x_coor = [0.0:1.0:3.0  3.75:0.25:6.25  7.0:1.0:10.0];
%y_coor = [0.0:1.0:2.0  2.25:0.25:3.75  4.0:1.0:6.0];
draw_mesh_alpha(x_coor,y_coor,z_coor,plotx,ploty,plotz,plotopt,gendata,savefilename);


