function initSolver_CNT()

display('Start solver initialization');

global kx ky kz scl;
global nodes links contacts;
global Nnode Nlink Nvolume; 
global nodeM volumeM;
global bndNodes intNodes dirNodes eqnNodes edgeNodes;
global eqnLinks  bndLinks ;
global isBndNodes isDirNodes dcVolDirNodes acVolDirNodes nodeLinks volumeNodes;
global doping;

global nodeLinks;
global linkL;

global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global sQMlinks;
global isqmnodes;
global currdlink;
global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;
global isqmvolm;
global metalNodes metalLinks;

%%%
global qkx qky qkz;
global dirxLinks diryLinks dirzLinks;
global JLinks J_amp;
global lightdirection Posz;
%%%

%%%%%%%%% define boundary and edge nodes %%%%%%%%%%%
[bndNodes,intNodes] = defBoundaryNodes(kx+1,ky+1,kz+1);	%bndNodes on 6 faces of the whole cube
isBndNodes = false(Nnode,1);
isBndNodes(bndNodes) = true;
edgeNodes = [];
for i = bndNodes'
    ajnd = nodeLinks{i}(2,:); %connecting nodes to bndNodes
    if all(isBndNodes(ajnd))
        edgeNodes = [edgeNodes;i]; %edge nodes on the 12 edge lines
    end
end


%%%%%%%%%% structure definition (in terms of node index) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%% plasmonic metal block %%%%%%%%%%%%%%%%%%
x_coor=[0:2.0:40];
y_coor=[0:2.0:40];
z_coor=[0:2.0:40];

Origin=[20.0,20.0,20.0];
[Xcen,Ycen,Zcen]=RealToNodes(Origin,x_coor,y_coor,z_coor);

% graphene is layed on x direction 
[X1,Y1,Z1]=RealToNodes([10.0,10.0,10.0],x_coor,y_coor,z_coor);
[X2,Y2,Z2]=RealToNodes([30.0,30.0,30.0],x_coor,y_coor,z_coor);

%[MX11,MY11,MZ11]=RealToNodes([22.0,26.0,26.0],x_coor,y_coor,z_coor);
%[MX12,MY12,MZ12]=RealToNodes([24.0,27.0,27.0],x_coor,y_coor,z_coor);

%[MX21,MY21,MZ21]=RealToNodes([29.0,26.0,26.0],x_coor,y_coor,z_coor);
%[MX22,MY22,MZ22]=RealToNodes([31.0,27.0,27.0],x_coor,y_coor,z_coor);

%display([' MX11 :',num2str(MX11),' MY11: ',num2str(MY11),' MZ11: ',num2str(MZ11)]);
%display([' MX12 :',num2str(MX12),' MY12: ',num2str(MY12),' MZ12: ',num2str(MZ12)]);
%display([' MX21 :',num2str(MX21),' MY21: ',num2str(MY21),' MZ21: ',num2str(MZ21)]);
%display([' MX22 :',num2str(MX22),' MY22: ',num2str(MY22),' MZ22: ',num2str(MZ22)]);
% x2-x1 = dev2ce x, y2-y1 = device y
[qmx1,qmy1,qmz1]=RealToNodes([10.0,10.0,10.0],x_coor,y_coor,z_coor);
[qmx2,qmy2,qmz2]=RealToNodes([30.0,30.0,30.0],x_coor,y_coor,z_coor);
%QM grid from lodestar
qkx =129 ; qky = 33; qkz = 33;
%qmx1 = 18; qmx2 = 25; qmy1 = 8; qmy2 = 22; qmz1 = 8; qmz2 = 22;
%R: set the boundary voltage to 0 now, for light
avoltage   = 0.0;
avoltageac = 0.0;
contacts{1,1} = defContacts([1,1,1],[1,ky+1,kz+1],avoltage/scl.Vt, avoltageac/scl.Vt);
contacts{2,1} = defContacts([kx+1,1,1],[kx+1,ky+1,kz+1],0,0); 
contacts{3,1} = defContacts([1,1,1],[kx+1,ky+1,1],0,0); 
contacts{4,1} = defContacts([1,1,kz+1],[kx+1,ky+1,kz+1],0,0); 
contacts{5,1} = defContacts([1,1,1],[kx+1,1,kz+1],0,0); 
contacts{6,1} = defContacts([1,ky+1,1],[kx+1,ky+1,kz+1],0,0);

dirNodes = [contacts{1}.nodes;contacts{2}.nodes;contacts{3}.nodes;contacts{4}.nodes;contacts{5}.nodes;contacts{6}.nodes];

isDirNodes = zeros(Nnode,1);
isDirNodes(dirNodes) = 1:length(dirNodes); %counter of dirNodes
dcVolDirNodes = [contacts{1}.dcVol;contacts{2}.dcVol;contacts{3}.dcVol;contacts{4}.dcVol;contacts{5}.dcVol;contacts{6}.dcVol];
acVolDirNodes = [contacts{1}.acVol;contacts{2}.acVol;contacts{3}.acVol;contacts{4}.acVol;contacts{5}.acVol;contacts{6}.acVol];
%Vector potential links
switch lightdirection
case 'kxEyBz'
	veclayer= defBrickLinks([kx+1,1,1],[kx+1,ky+1,kz+1]);
	JLinks=intersect(veclayer,diryLinks);
case 'kxEzBy'
	veclayer= defBrickLinks([kx+1,1,1],[kx+1,ky+1,kz+1]);
	JLinks=intersect(veclayer,dirzLinks);
case 'kyExBz'
	veclayer= defBrickLinks([1,ky+1,1],[kx+1,ky+1,kz+1]);
	JLinks=intersect(veclayer,dirxLinks);
case 'kyEzBx'
	veclayer= defBrickLinks([1,ky+1,1],[kx+1,ky+1,kz+1]);
	JLinks=intersect(veclayer,dirzLinks);
case 'kzEyBx'	
	veclayer= defBrickLinks([1,1,kz+1],[kx+1,ky+1,kz+1]);
	JLinks=intersect(veclayer,diryLinks);
case 'kzExBy'	
	veclayer= defBrickLinks([1,1,kz+1],[kx+1,ky+1,kz+1]);
	JLinks=intersect(veclayer,dirxLinks);
case 'nolight'
	JLinks=[];
case 'specialkzExBy'
	tmpz=zRealTozNodes(Posz,z_coor);
	veclayer= defBrickLinks([1,1,tmpz],[kx+1,ky+1,tmpz]);
	JLinks=intersect(veclayer,dirxLinks);
case 'specialkzEyBx'	
	tmpz=zRealTozNodes(Posz,z_coor);
	veclayer= defBrickLinks([1,1,tmpz],[kx+1,ky+1,tmpz]);
	JLinks=intersect(veclayer,diryLinks);
end
J_amp=6.0*1e6/(40*1e-9)*ones(length(JLinks),1)/scl.s_J;

%%%%%%%
Radius=10.0;
metalNodes=defSphereNodes([20.0,20.0,20.0],scl.lambda,Radius,nodes);
%%%%%%%
savefilename='metalNodes.mat';
save(savefilename, 'metalNodes','Nnode','kx','ky','kz','x_coor','y_coor','z_coor');

%QMnodes = [];
%sQMnodes = [];
 QMnodes =defBrickNodes([qmx1,qmy1,qmz1],[qmx2,qmy2,qmz2]);
 sQMnodes =defBrickNodes([qmx1+1,qmy1+1,qmz1+1],[qmx2-1,qmy2-1,qmz2-1]);

istqmnodes            = false(Nnode,1);
istqmnodes(QMnodes)   = true;

isSqmnodes            = false(Nnode,1);
isSqmnodes(sQMnodes)  = true;
isqmnodes             = false(Nnode,1);
isqmnodes(sQMnodes)   = true;
isSQMlinks            = zeros(Nlink,1);

QMlinks   =[];
sQMlinks  =[];
currdlink = zeros(Nlink,6);


for ilk =1:Nlink

    n1  = links(ilk,1); n2 = links(ilk,2); n3 = links(ilk,3);

    [nx,ny,nz]=id2co(n1); [nx1,ny1,nz1]=id2co(n2);

    currdlink(ilk,1)=(nx+nx1)/2;
    currdlink(ilk,2)=(ny+ny1)/2;
    currdlink(ilk,3)=(nz+nz1)/2; 

    if (isSqmnodes(n1)+isSqmnodes(n2))==1
      sQMlinks = [sQMlinks;ilk];
       QMlinks = [QMlinks;ilk];
      if n3==1
         if nx < ( X1 + X2 ) /2
          isSQMlinks(ilk)  = 1;
         else
          isSQMlinks(ilk)  = 2;
         end
      elseif n3==2 
         if ny < ( Y1 + Y2 ) /2
          isSQMlinks(ilk)  = 3;        
         else
          isSQMlinks(ilk)  = 4;
         end
      elseif n3==3 
         if nz < ( Z1 + Z2 ) /2
          isSQMlinks(ilk)  = 5;         
         else 
          isSQMlinks(ilk)  = 6;
         end
      end 
    end

    if (isqmnodes(n1)+isqmnodes(n2))==2
      QMlinks = [QMlinks;ilk];
      if n3==1
          isSQMlinks(ilk)  = 7;
      elseif n3==2
          isSQMlinks(ilk)  = 8;
      elseif n3==3
          isSQMlinks(ilk)  = 9;
      end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idMetalNodes = false(Nnode,1);
idMetalNodes(metalNodes) = true;
%idInsNodes = false(Nnode,1);
%idInsNodes(insNodes) = true;

%%%% define doping profile %%%%%%%%%%%%%%

%%%%%%%%%%%%
bndLinks = [];
for i = 1:Nlink
    n1 = links(i,1); n2 = links(i,2);
    if isBndNodes(n1) && isBndNodes(n2)
       bndLinks = [bndLinks;i]; 
    end
end

metalLinks=[];
for i = 1:Nlink
    n1 = links(i,1); n2 = links(i,2);
    if idMetalNodes(n1) && idMetalNodes(n2)
       metalLinks = [metalLinks;i]; 
    end
end

%%% material types of each volume (semiconductor = 1, metal = 2, insulator =3 %%%
%%% second insulator  = 4		
volumeM = 2*ones(Nvolume,1);
for i = 1:Nvolume
   vNodes = volumeNodes(i,:);
   if all(idMetalNodes(vNodes))
      volumeM(i) = 1; 
   end
end

isqmvolm = false(Nvolume,1);
for i = 1:Nvolume
   vNodes = volumeNodes(i,:);
   if all(istqmnodes(vNodes))
      isqmvolm(i) = true;
   end
end

display('End solver initialization');