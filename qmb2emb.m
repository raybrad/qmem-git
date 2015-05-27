function  qmb2emb()

global kx ky kz;
global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;
%
global qkx qky qkz;
global currdfx currdfy currdfz;
global currdemx currdemy currdemz;
global currdx currdy currdz;

ix = qmx2-qmx1+1; iy = qmy2-qmy1+1; iz = qmz2-qmz1+1;

qmx = 4; qmy = 1; qmz=1;

load('qmboundcd');

qkx = qmboundcd(1,1); qky = qmboundcd(1,2); qkz = qmboundcd(1,3); 

qmboundcd(1,:)=[];

currdx = reshape(qmboundcd(:,1),qkx,qky,qkz); 
currdy = reshape(qmboundcd(:,2),qkx,qky,qkz);
currdz = reshape(qmboundcd(:,3),qkx,qky,qkz);

[xi,yi,zi] = meshgrid(0:qmy/(qky-1):qmy,0:qmx/(qkx-1):qmx,0:qmz/(qkz-1):qmz);
[x,y,z]=meshgrid(0:qmy/(iy-1):qmy,0:qmx/(ix-1):qmx,0:qmz/(iz-1):qmz);

currdemx=interp3(xi,yi,zi,currdx,x,y,z,'spline');  
currdemy=interp3(xi,yi,zi,currdy,x,y,z,'spline');
currdemz=interp3(xi,yi,zi,currdz,x,y,z,'spline');

currdfx=reshape(currdemx,ix*iy*iz,1);
currdfy=reshape(currdemy,ix*iy*iz,1);
currdfz=reshape(currdemz,ix*iy*iz,1);








