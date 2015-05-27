function draw_mesh_alpha(x_coor,y_coor,z_coor,plotx,ploty,plotz,plotopt,gendata,savefilename)
% A simple code to plot and generate 3D, cartesian, non voronoi, mesh
%{
%  input: (see examples below)
%        x_coor      : nodes on x-axis
%        y_coor      : nodes on y-axis
%        z_coor      : nodes on z-axis
%        plotx       : plot the grids on surfaces normal to x-axis
%        ploty       : plot the grids on surfaces normal to y-axis
%        plotz       : plot the grids on surfaces normal to z-axis
%        plotopt     : mode of plotting
%        gendata     : generate and output mesh data
%        savefilename: name of the file saving the mesh data
%}
%
%clc;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

%input examples:
%x_coor = [0.0:1.0:3.0  3.75:0.25:6.25  7.0:1.0:10.0];
%y_coor = [0.0:1.0:2.0  2.25:0.25:3.75  4.0:1.0:6.0];
%z_coor = [0  0.75:0.25:2.25  8/3  10/3  3.75:0.25:5.25  6];


%plotx = true;
%ploty = true;
%plotz = true;

%plotopt=2;  %1: plotting grid for the whole volume; 
            %2: plotting the first surface along each direction
            %3: plotting two end surfaces along each direction

%gendata = true;  % true  = generate and output mesh data
                 % false = NOT do the work

%savefilename = 'test_mesh.mat';   %the output file will be in .MAT format

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x_coor = [0.0:1.0:3.0  3.75:0.25:6.25  7.0:1.0:10.0];
%y_coor = [0.0:1.0:2.0  2.25:0.25:3.75  4.0:1.0:6.0];
%z_coor = [0  0.75:0.25:2.25  8/3  10/3  3.75:0.25:5.25  6];

%x_coor = [0.0:0.5:6.0];
%y_coor = [0.0:0.5:12.0];
%z_coor = [0.0:0.5:6.0];

no_of_nodes_x = length(x_coor);   %number of nodes along x-axis
no_of_nodes_y = length(y_coor);   %number of nodes along y-axis
no_of_nodes_z = length(z_coor);   %number of nodes along z-axis

Nx=no_of_nodes_x;
Ny=no_of_nodes_y;
Nz=no_of_nodes_z;

Nnode=Nx*Ny*Nz;
Nlink=3*Nx*Ny*Nz-Nx*Ny-Nx*Nz-Ny*Nz;
Nsurf=3*Nx*Ny*Nz-2*Nx*Ny-2*Nx*Nz-2*Ny*Nz+Nx+Ny+Nz;
Nvolume=(Nx-1)*(Ny-1)*(Nz-1);

display([' Nx: ',num2str(no_of_nodes_x),' Ny: ',num2str(no_of_nodes_y),' Nz :',num2str(no_of_nodes_z)]);
display([' Nnode: ',num2str(Nnode)]);	
display([' Nlink: ',num2str(Nlink)]);	
display([' Nsurf: ',num2str(Nsurf)]);	
display([' Nvolume: ',num2str(Nvolume)]);	

drawgrid=zeros(4*(Nx*Ny+Nx*Nz+Ny*Nz),3);
nodes   =zeros(Nnode,4);
tmpnodes=zeros(Nx,Ny,Nz);
links   =zeros(Nlink,3);
templink=zeros(2*Nx,2*Ny,2*Nz);
surfLinks =zeros(Nsurf*4,2);
tempsurf  =zeros(2*Nx,2*Ny,2*Nz);
volumeNodes=zeros(Nvolume*8,3);
volumeLinks=zeros(Nvolume*12,3);
volumeSurfs=zeros(Nvolume*6,3);	
%nodeV(Nnode),linkSs(Nlink),dlinkL(Nsurf),nodeVolV(Nnode,Nvolume),linkVolS(Nlink,Nvolume))

if ~(plotx || ploty || plotz)
    display('Warning! No grids are plotting!');
%    return
else
    tmpcounter = 0;
    switch plotopt
        case 1
            x_loop = 1:no_of_nodes_x;
            y_loop = 1:no_of_nodes_y;
            z_loop = 1:no_of_nodes_z;
        case 2
            x_loop = 1:1;
            y_loop = 1:1;
            z_loop = 1:1;
        case 3
            x_loop = 1:no_of_nodes_x-1:no_of_nodes_x;
            y_loop = 1:no_of_nodes_y-1:no_of_nodes_y;
            z_loop = 1:no_of_nodes_z-1:no_of_nodes_z;
        otherwise
            display('Error! The value of plotopt should be 1, 2 or 3.');
    end
end

% along x-direction
if plotx
    for i = x_loop
        for j = 1:no_of_nodes_y
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(i);
            drawgrid(tmpcounter,2)= y_coor(j);
            drawgrid(tmpcounter,3)= z_coor(1);
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(i);
            drawgrid(tmpcounter,2)= y_coor(j);
            drawgrid(tmpcounter,3)= z_coor(no_of_nodes_z);
        end
        for k = 1:no_of_nodes_z
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(i);
            drawgrid(tmpcounter,2)= y_coor(1);
            drawgrid(tmpcounter,3)= z_coor(k);
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(i);
            drawgrid(tmpcounter,2)= y_coor(no_of_nodes_y);
            drawgrid(tmpcounter,3)= z_coor(k);
        end
    end
end
% along y-direction
if ploty
    for j = y_loop
        for i = 1:no_of_nodes_x
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(i);
            drawgrid(tmpcounter,2)= y_coor(j);
            drawgrid(tmpcounter,3)= z_coor(1);
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(i);
            drawgrid(tmpcounter,2)= y_coor(j);
            drawgrid(tmpcounter,3)= z_coor(no_of_nodes_z);
        end
        for k = 1:no_of_nodes_z
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(1);
            drawgrid(tmpcounter,2)= y_coor(j);
            drawgrid(tmpcounter,3)= z_coor(k);
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(no_of_nodes_x);
            drawgrid(tmpcounter,2)= y_coor(j);
            drawgrid(tmpcounter,3)= z_coor(k);
        end
    end
end
% along z-direction
if plotz
    for k = z_loop
        for j = 1:no_of_nodes_y
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(1);
            drawgrid(tmpcounter,2)= y_coor(j);
            drawgrid(tmpcounter,3)= z_coor(k);
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(no_of_nodes_x);
            drawgrid(tmpcounter,2)= y_coor(j);
            drawgrid(tmpcounter,3)= z_coor(k);
        end
        for i = 1:no_of_nodes_x
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(i);
            drawgrid(tmpcounter,2)= y_coor(1);
            drawgrid(tmpcounter,3)= z_coor(k);
            tmpcounter = 1 + tmpcounter;
            drawgrid(tmpcounter,1)= x_coor(i);
            drawgrid(tmpcounter,2)= y_coor(no_of_nodes_y);
            drawgrid(tmpcounter,3)= z_coor(k);
        end
    end
end

%plotting the grids
%figure;
%hold on;
%xlabel('X - axis');
%ylabel('Y - axis');
%zlabel('Z - axis');
%for i = 1:2:tmpcounter
%    plot3(drawgrid(i:i+1,1),drawgrid(i:i+1,2),drawgrid(i:i+1,3))
%end

%{--------------------Mesh Data Generation Starts Here------------------ %}
if (gendata)
%
    display('Gen data begin');
    % Defining nodes
    % nodes = nodes(order number, x-coor, y-coor, z-coor)
    tmpcounter = 0;
    for k = 1:no_of_nodes_z
        for j = 1:no_of_nodes_y
            for i = 1:no_of_nodes_x
    %            tmpcounter = (no_of_nodes_x*(j-1))+(no_of_nodes_y*no_of_nodes_x*(k-1))+i;
                tmpcounter = 1 + tmpcounter;
                nodes(tmpcounter,1) = tmpcounter;
                nodes(tmpcounter,2) = x_coor(i);
                nodes(tmpcounter,3) = y_coor(j);
                nodes(tmpcounter,4) = z_coor(k);
                tmpnodes(i,j,k) = tmpcounter;
            end
        end
    end
    %
    display('Gen nodes');

    % Defining links
    % Two nodes form a link
    % links = links(order number, starting node, ending node)
    tmpcounter = 0;    
    for k = 1:no_of_nodes_z
        for j = 1:no_of_nodes_y
            for i = 1:no_of_nodes_x
                if (i < no_of_nodes_x) 
                tmpcounter = 1 + tmpcounter;
                links(tmpcounter,1) = tmpcounter;
                links(tmpcounter,2) = tmpnodes(i,j,k);
                links(tmpcounter,3) = tmpnodes(i+1,j,k);
                templink((2*i),(2*j)-1,(2*k)-1) = tmpcounter;
                end
                if (j < no_of_nodes_y) 
                tmpcounter = 1 + tmpcounter;
                links(tmpcounter,1) = tmpcounter;
                links(tmpcounter,2) = tmpnodes(i,j,k);
                links(tmpcounter,3) = tmpnodes(i,j+1,k);
                templink((2*i)-1,(2*j),(2*k)-1) = tmpcounter;
                end
                if (k < no_of_nodes_z) 
                tmpcounter = 1 + tmpcounter;
                links(tmpcounter,1) = tmpcounter;
                links(tmpcounter,2) = tmpnodes(i,j,k);
                links(tmpcounter,3) = tmpnodes(i,j,k+1);
                templink((2*i)-1,(2*j)-1,(2*k)) = tmpcounter;
                end
            end
        end
    end
    %
    display('Gen links');

    % Defining surfaces with links
    % four links form a surface
    % surfLinks = surfLinks(order number, one of the four surrounding links)
    tmpcounter = 0;
    tmpcounter1 = 0;
    for k = 1:no_of_nodes_z
        for j = 1:no_of_nodes_y
            for i = 1:no_of_nodes_x
                if ((i < no_of_nodes_x) && (j < no_of_nodes_y))
                    tmpcounter = 1 + tmpcounter;
                    tmpcounter1 = 1 + tmpcounter1;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i),(2*j)-1,(2*k)-1);
                    tmpcounter = 1 + tmpcounter;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i)+1,(2*j),(2*k)-1);
                    tmpcounter = 1 + tmpcounter;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i),(2*j)+1,(2*k)-1);
                    tmpcounter = 1 + tmpcounter;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i)-1,(2*j),(2*k)-1);
                    tempsurf((2*i),(2*j),(2*k)-1) = tmpcounter1;
                end
                if ((j < no_of_nodes_y) && (k < no_of_nodes_z))
                    tmpcounter = 1 + tmpcounter;
                    tmpcounter1 = 1 + tmpcounter1;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i)-1,(2*j),(2*k)-1);
                    tmpcounter = 1 + tmpcounter;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i)-1,(2*j)+1,(2*k));
                    tmpcounter = 1 + tmpcounter;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i)-1,(2*j),(2*k)+1);
                    tmpcounter = 1 + tmpcounter;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i)-1,(2*j)-1,(2*k));
                    tempsurf((2*i)-1,(2*j),(2*k)) = tmpcounter1;
                end
                if ((k < no_of_nodes_z) && (i < no_of_nodes_x))
                    tmpcounter = 1 + tmpcounter;
                    tmpcounter1 = 1 + tmpcounter1;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i),(2*j)-1,(2*k)-1);
                    tmpcounter = 1 + tmpcounter;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i)+1,(2*j)-1,(2*k));
                    tmpcounter = 1 + tmpcounter;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i),(2*j)-1,(2*k)+1);
                    tmpcounter = 1 + tmpcounter;
                    surfLinks(tmpcounter, 1) = tmpcounter1;
                    surfLinks(tmpcounter, 2) = templink((2*i)-1,(2*j)-1,(2*k));
                    tempsurf((2*i),(2*j)-1,(2*k)) = tmpcounter1;
                end
            end
        end
    end
    %
    display('Gen surflinks');

    % Defining volume with nodes
    % eight nodes form a volume
    % volumeNodes = volumeNodes(order number, one of the eight surrounding
    % nodes, one eighth of the volume of the corresponding cube)
    tmpcounter = 0;
    tmpcounter1 = 0;
    for k = 1:(no_of_nodes_z-1)
        for j = 1:(no_of_nodes_y-1)
            for i = 1:(no_of_nodes_x-1)
                tmp =  (nodes(tmpnodes(i+1,j,k),2) - nodes(tmpnodes(i,j,k),2)) ...
                        * (nodes(tmpnodes(i,j+1,k),3) - nodes(tmpnodes(i,j,k),3)) ...
                        * (nodes(tmpnodes(i,j,k+1),4) - nodes(tmpnodes(i,j,k),4)) / 8;
                tmpcounter = 1 + tmpcounter;
                tmpcounter1 = 1 + tmpcounter1;
                volumeNodes(tmpcounter, 1) = tmpcounter1;
                volumeNodes(tmpcounter, 2) = tmpnodes(i,j,k);
                volumeNodes(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeNodes(tmpcounter, 1) = tmpcounter1;
                volumeNodes(tmpcounter, 2) = tmpnodes(i+1,j,k);
                volumeNodes(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeNodes(tmpcounter, 1) = tmpcounter1;
                volumeNodes(tmpcounter, 2) = tmpnodes(i,j+1,k);
                volumeNodes(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeNodes(tmpcounter, 1) = tmpcounter1;
                volumeNodes(tmpcounter, 2) = tmpnodes(i,j,k+1);
                volumeNodes(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeNodes(tmpcounter, 1) = tmpcounter1;
                volumeNodes(tmpcounter, 2) = tmpnodes(i+1,j+1,k);
                volumeNodes(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeNodes(tmpcounter, 1) = tmpcounter1;
                volumeNodes(tmpcounter, 2) = tmpnodes(i+1,j,k+1);
                volumeNodes(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeNodes(tmpcounter, 1) = tmpcounter1;
                volumeNodes(tmpcounter, 2) = tmpnodes(i,j+1,k+1);
                volumeNodes(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeNodes(tmpcounter, 1) = tmpcounter1;
                volumeNodes(tmpcounter, 2) = tmpnodes(i+1,j+1,k+1);
                volumeNodes(tmpcounter, 3) = tmp;
            end
        end
    end
    %
    display('Gen volumenodes');

    % Defining volume with links
    % twelve links form a volume
    % volumeLinks = volumeLinks(order number, one of the twelve surrounding
    % links, one fourth of the area of the orthogonal surface in the
    % corresponding cube)
    tmpcounter = 0;
    tmpcounter1 = 0;
    for k = 1:(no_of_nodes_z-1)
        for j = 1:(no_of_nodes_y-1)
            for i = 1:(no_of_nodes_x-1)
                tmpcounter1 = 1 + tmpcounter1;
                tmp =  (nodes(tmpnodes(i,j+1,k),3) - nodes(tmpnodes(i,j,k),3)) ...
                        * (nodes(tmpnodes(i,j,k+1),4) - nodes(tmpnodes(i,j,k),4)) / 4;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i),(2*j)-1,(2*k)-1);
                volumeLinks(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i),(2*j)-1,(2*k)+1);
                volumeLinks(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i),(2*j)+1,(2*k)-1);
                volumeLinks(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i),(2*j)+1,(2*k)+1);
                volumeLinks(tmpcounter, 3) = tmp;
                tmp =  (nodes(tmpnodes(i+1,j,k),2) - nodes(tmpnodes(i,j,k),2)) ...
                        * (nodes(tmpnodes(i,j,k+1),4) - nodes(tmpnodes(i,j,k),4)) / 4;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i)-1,(2*j),(2*k)-1);
                volumeLinks(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i)-1,(2*j),(2*k)+1);
                volumeLinks(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i)+1,(2*j),(2*k)-1);
                volumeLinks(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i)+1,(2*j),(2*k)+1);
                volumeLinks(tmpcounter, 3) = tmp;
                tmp =  (nodes(tmpnodes(i,j+1,k),3) - nodes(tmpnodes(i,j,k),3)) ...
                        * (nodes(tmpnodes(i+1,j,k),2) - nodes(tmpnodes(i,j,k),2)) / 4;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i)-1,(2*j)-1,(2*k));
                volumeLinks(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i)-1,(2*j)+1,(2*k));
                volumeLinks(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i)+1,(2*j)-1,(2*k));
                volumeLinks(tmpcounter, 3) = tmp;
                tmpcounter = 1 + tmpcounter;
                volumeLinks(tmpcounter, 1) = tmpcounter1;
                volumeLinks(tmpcounter, 2) = templink((2*i)+1,(2*j)+1,(2*k));
                volumeLinks(tmpcounter, 3) = tmp;
            end
        end
    end
    %
    display('Gen volumelinks');

    % Defining volume with surfaces
    % six surfaces form a volume
    % volumeSurfs = volumeSurfs(order number, one of the six surrounding
    % surfaces, one half of the prependicular length in the corresponding cube)
    tmpcounter = 0;
    tmpcounter1 = 0;
    for k = 1:(no_of_nodes_z-1)
        for j = 1:(no_of_nodes_y-1)
            for i = 1:(no_of_nodes_x-1)
                tmpcounter = 1 + tmpcounter;
                tmpcounter1 = 1 + tmpcounter1;
                volumeSurfs(tmpcounter, 1) = tmpcounter1;
                volumeSurfs(tmpcounter, 2) = tempsurf((2*i),(2*j),(2*k)-1);
                volumeSurfs(tmpcounter, 3) = (nodes(tmpnodes(i,j,k+1),4) - nodes(tmpnodes(i,j,k),4)) / 2;
                tmpcounter = 1 + tmpcounter;
                volumeSurfs(tmpcounter, 1) = tmpcounter1;
                volumeSurfs(tmpcounter, 2) = tempsurf((2*i)-1,(2*j),(2*k));
                volumeSurfs(tmpcounter, 3) = (nodes(tmpnodes(i+1,j,k),2) - nodes(tmpnodes(i,j,k),2)) / 2;
                tmpcounter = 1 + tmpcounter;
                volumeSurfs(tmpcounter, 1) = tmpcounter1;
                volumeSurfs(tmpcounter, 2) = tempsurf((2*i),(2*j)-1,(2*k));
                volumeSurfs(tmpcounter, 3) = (nodes(tmpnodes(i,j+1,k),3) - nodes(tmpnodes(i,j,k),3)) / 2;
                tmpcounter = 1 + tmpcounter;
                volumeSurfs(tmpcounter, 1) = tmpcounter1;
                volumeSurfs(tmpcounter, 2) = tempsurf((2*i)+1,(2*j),(2*k));
                volumeSurfs(tmpcounter, 3) = (nodes(tmpnodes(i+1,j,k),2) - nodes(tmpnodes(i,j,k),2)) / 2;
                tmpcounter = 1 + tmpcounter;
                volumeSurfs(tmpcounter, 1) = tmpcounter1;
                volumeSurfs(tmpcounter, 2) = tempsurf((2*i),(2*j)+1,(2*k));
                volumeSurfs(tmpcounter, 3) = (nodes(tmpnodes(i,j+1,k),3) - nodes(tmpnodes(i,j,k),3)) / 2;
                tmpcounter = 1 + tmpcounter;
                volumeSurfs(tmpcounter, 1) = tmpcounter1;
                volumeSurfs(tmpcounter, 2) = tempsurf((2*i),(2*j),(2*k)+1);
                volumeSurfs(tmpcounter, 3) = (nodes(tmpnodes(i,j,k+1),4) - nodes(tmpnodes(i,j,k),4)) / 2;
            end
        end
    end
    %
    display('gen volumesurfs');
    save(savefilename, 'links', 'nodes', 'surfLinks', 'volumeLinks', 'volumeNodes', 'volumeSurfs');
    %
end
    display('Good! Job Completed!');
