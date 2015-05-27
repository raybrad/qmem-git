function plotStructureIso
%load 'metalNodes.mat';
load 'metalNodes2.mat';

Nnode=250047;
kx=62;
ky=62;
kz=62;
x_coor=[0:1:9  9.5:0.5:30.5 31:1:40];
y_coor=[0:1:9  9.5:0.5:30.5 31:1:40];
z_coor=[0:1:9  9.5:0.5:30.5 31:1:40];
ismetalNodes=zeros(Nnode,1);
ismetalNodes(metalNodes)=1;

nodex=kx+1;
nodey=ky+1;
nodez=kz+1;

isTranmetalNodes(1:nodex,1:nodey,1:nodez)=0;
tmpcounter=0;
for k = 1:nodez
    for j = 1:nodey
        for i = 1:nodex
	tmpcounter=tmpcounter+1;
	isTranmetalNodes(i,j,k)=ismetalNodes(tmpcounter);
	end
    end
end
[X,Y]=ndgrid(x_coor,y_coor);
%[X,Y]=ndgrid(1:1:68,1:1:68);
%aveE(:,:)=isTranmetalNodes(:,46,:);
aveE(:,:)=isTranmetalNodes(:,:,32);
figure;
contourf(X,Y,aveE,1);
%contour(X,Y,aveE,20);

%scatter3(X,Y,Z,isTranmetalNodes);
%exit;
clear;
