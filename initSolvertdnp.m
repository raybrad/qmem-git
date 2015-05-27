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

global semiNodes eqnSemiNodes;
global dirSemiNodes isDirSemiNodes;
global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global sQMlinks;
global isqmnodes;
global currdlink;
global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;
global isqmvolm;
global issemvolm_p;
global semiNodesp semiNodesn;
global semiintf;
global issemiintf;
global npvaryped;
global seminodevolm seminodevolmn seminodevolmp;

%%%%%%%%% define boundary and edge nodes %%%%%%%%%%%
[bndNodes,intNodes] = defBoundaryNodes(kx+1,ky+1,kz+1);
isBndNodes = false(Nnode,1);
isBndNodes(bndNodes) = true;
edgeNodes = [];
for i = bndNodes'
    ajnd = nodeLinks{i}(2,:);
    if all(isBndNodes(ajnd))
        edgeNodes = [edgeNodes;i];
    end
end


%%%%%%%%%% structure definition (in terms of node index) %%%%%%%%%%%%%%%%%%
%%%%%%%%% CNT0 %%%%%%%%%%%%%%%%%
% X1 = 4; X2 = 14; Y1 = 4; Y2 = 14; Z1 = 4; Z2 = 14;
% 
% contacts{1,1} = defContacts([1,Y1,Z1],[1,Y2,Z2],0.1/scl.Vt,0.1/scl.Vt);
% contacts{2,1} = defContacts([kx+1,Y1,Z1],[kx+1,Y2,Z2],0,0);
% dirNodes = [contacts{1}.nodes;contacts{2}.nodes];
% isDirNodes = zeros(Nnode,1);
% isDirNodes(dirNodes) = 1:length(dirNodes);
% dcVolDirNodes = [contacts{1}.dcVol;contacts{2}.dcVol];
% acVolDirNodes = [contacts{1}.acVol;contacts{2}.acVol];
% 
% semiNodes = defBrickNodes([X1,Y1,Z1],[X2,Y2,Z2]); % semiconductor region
% metalNodes1 = defBrickNodes([1,Y1,Z1],[X1,Y2,Z2]);
% metalNodes2 = defBrickNodes([X2,Y1,Z1],[kx+1,Y2,Z2]);
% metalNodes = [metalNodes1;metalNodes2]; % metal region
% QMnodes = defBrickNodes([X1+2,Y1+2,Z1+2],[X2-2,Y2-2,Z2-2]);
% 
% interfNodes = intersect(semiNodes,metalNodes); % interface between metal and semiconductor
% dirSemiNodes = interfNodes;
% isDirSemiNodes = false(Nnode,1);
% isDirSemiNodes(dirSemiNodes) = true;

%%%%%%%%%%% CNT1 %%%%%%%%%%%%%%%%%%
X1 = 4; X2 = 44; Y1 = 3; Y2 = 19; Z1 = 3; Z2 = 19; Xm1 = 17; Xm2 = 13; Ym = 9; Ym2 = 21; Zm = 9;
Ym1 = 9;
qmx1 = X1; qmx2 = X2; qmy1 = Y1; qmy2 = Y2; qmz1 = Z1; qmz2 = Z2;
avoltage1  = 1.0e-0;
avoltage2  = 0.0e-3;
%avoltageg  = -1.0e-3;

contacts{1,1} = defContacts([1,Y1,Z1],[1,Y2,Z2],avoltage1/scl.Vt,0);
contacts{2,1} = defContacts([kx+1,Y1,Z1],[kx+1,Y2,Z2],avoltage2/scl.Vt,0);
%contacts{3,1} = defContacts([Xm1+1,2,kz+1],[Xm2-1,ky,kz+1],0,avoltageg/scl.Vt);
%contacts{5,1} = defContacts([Xm1+1,ky+1,Z1+1],[Xm2-1,ky+1,Z2-1],avoltageg/scl.Vt,0);
%contacts{3,1} = defContacts([X1+1,Y1+1,kz+1],[X2-1,Y2-1,kz+1],0,0);
dirNodes = [contacts{1}.nodes;contacts{2}.nodes];
%dirNodes = [contacts{1}.nodes;contacts{2}.nodes;contacts{3}.nodes];
isDirNodes = zeros(Nnode,1);
isDirNodes(dirNodes) = 1:length(dirNodes);
dcVolDirNodes = [contacts{1}.dcVol;contacts{2}.dcVol];
acVolDirNodes = [contacts{1}.acVol;contacts{2}.acVol];
%dcVolDirNodes = [contacts{1}.dcVol;contacts{2}.dcVol;contacts{3}.dcVol];
%acVolDirNodes = [contacts{1}.acVol;contacts{2}.acVol;contacts{3}.acVol];

%semiNodes1 = defBrickNodes([1,Y1,Z1],[kx+1,Y2,Z2]);
%semiNodes2 = defBrickNodes([X1,Y1,Z1],[X2,Y2,kz+1]);
%semiNodes3 = defBrickNodes([X1,Y1,Z1],[X2,Y2,Z2]); % semiconductor region

%semiNodes1 = defBrickNodes([1,1,1],[kx+1,ky+1,Z1]);
%semiNodes2 = defBrickNodes([1,1,Z1],[kx+1,Y1,Z2]); 
%semiNodes3 = defBrickNodes([X1,Y1,Z1],[X2,Y2,Z2]); % semiconductor region
%semiNodes4 = defBrickNodes([1,Y2,Z1],[kx+1,ky+1,Z2]);
%semiNodes5 = defBrickNodes([1,1,Z2],[kx+1,ky+1,kz+1]);
semiNodes1=[];
semiNodes2=[];

seminodevolmn=[];
seminodevolmp=[];

for ddz = Z1:2*npvaryped:Z2-2*npvaryped
  for ddy = Y1:2*npvaryped:Y2-2*npvaryped

   semiNodes_nnp = defBrickNodes([X1,ddy,ddz],[X2,(ddy+npvaryped),(ddz+npvaryped)]);
   semiNodes1 = unique([semiNodes1;semiNodes_nnp]);
   semiNodes_vol = defBrickVolm([X1,ddy,ddz],[X2,(ddy+npvaryped),(ddz+npvaryped)]);
   seminodevolmp = unique([seminodevolmp;semiNodes_vol]);
  
   semiNodes_nnp = defBrickNodes([X1,(ddy+npvaryped),ddz],[X2,(ddy+2*npvaryped),(ddz+npvaryped)]);
   semiNodes2 = unique([semiNodes2;semiNodes_nnp]);
   semiNodes_vol = defBrickVolm([X1,(ddy+npvaryped),ddz],[X2,(ddy+2*npvaryped),(ddz+npvaryped)]);
   seminodevolmn = unique([seminodevolmn;semiNodes_vol]);

  end
end

for ddz = Z1+npvaryped : 2*npvaryped : Z2-npvaryped
  for ddy = Y1 : 2*npvaryped : Y2-2*npvaryped

   semiNodes_nnp = defBrickNodes([X1,ddy,ddz],[X2,(ddy+npvaryped),(ddz+npvaryped)]);
   semiNodes2 = unique([semiNodes2;semiNodes_nnp]);
   semiNodes_vol = defBrickVolm([X1,ddy,ddz],[X2,(ddy+npvaryped),(ddz+npvaryped)]);
   seminodevolmn = unique([seminodevolmn;semiNodes_vol]);

   semiNodes_nnp = defBrickNodes([X1,(ddy+npvaryped),ddz],[X2,(ddy+2*npvaryped),(ddz+npvaryped)]);
   semiNodes1 = unique([semiNodes1;semiNodes_nnp]);
   semiNodes_vol = defBrickVolm([X1,(ddy+npvaryped),ddz],[X2,(ddy+2*npvaryped),(ddz+npvaryped)]);
   seminodevolmp = unique([seminodevolmp;semiNodes_vol]);

  end
end

%semiNodes1 = defBrickNodes([X1,Y1,Z1],[X2,Ym,Zm]);
%semiNodes2 = defBrickNodes([X1,Ym,Zm],[X2,Y2,Z2]);

semiNodes3 = defBrickNodes([X1,Ym,Z1],[X2,Y2,Zm]);
%semiNodes4 = defBrickNodes([X1,Y1,Zm],[X2,Ym,Z2]);
semiNodes5 = defBrickNodes([1,Y1,Z1],[X1,Y2,Z2]);
semiNodes6 = defBrickNodes([X2,Y1,Z1],[kx+1,Y2,Z2]);
semiNodes_vol5 = defBrickVolm([1,Y1,Z1],[X1,Y2,Z2]);
semiNodes_vol6 = defBrickVolm([X2,Y1,Z1],[kx+1,Y2,Z2]);

%semiNodes6 = defBrickNodes([1,Y1,Z1],[X1,Y2,Z2]);
%semiNodes7 = defBrickNodes([X2,Y1,Z1],[kx+1,Y2,Z2]);

semiNodes10= defBrickNodes([X1+1,Y1+1,Z1+1],[X2-1,Y2-1,Z2-1]);

%semiNodes13= defBrickNodes([X1,Y1+1,Z1+1],[X2,Y2-1,Z2-1]);

semiNodes14 = defBrickNodes([1,Y1,Z1],[kx+1,Y2,Z2]);

semiNodesp = unique([semiNodes5;semiNodes1]);
semiNodesn = unique([semiNodes6;semiNodes2]);
seminodevolmp = unique([seminodevolmp;semiNodes_vol5]);
seminodevolmn = unique([seminodevolmn;semiNodes_vol6]);
semiNodes = unique([semiNodes14]);

issemcond_p = false(Nnode,1);
issemcond_p(semiNodesp) = true;

issemcond_n = false(Nnode,1);
issemcond_n(semiNodesn) = true;
%semiNodes = unique([semiNodes6;semiNodes7]);
%semiNodes = unique([semiNodes1;semiNodes2]);
%semiNodes = unique([semiNodes3]);
%semiNodes   = semiNodes3;
%semiNodes   = [];
%metalNodes1 = defBrickNodes([1,Y1,Z1],[X1,Ym1,Z2]);
%metalNodes2 = defBrickNodes([1,Ym2,Z1],[X1,Y2,Z2]);
%metalNodes3 = defBrickNodes([X2,Y1,Z1],[kx+1,Y2,Z2]);
%metalNodes4 = defBrickNodes([Xm1,Ym3,Z1],[Xm2,Ym4,kz+1]);
%metalNodes5 = setdiff(metalNodes4,[semiNodes7;semiNodes8]);
%metalNodes4 = defBrickNodes([X2,Y1,Zt1],[kx+1,Y2,Zt2]);
%metalNodes = [metalNodes1;metalNodes2]; % metal region
%metalNodes = [metalNodes1;metalNodes2;metalNodes3]; % metal regio
%metalNodes = unique([metalNodes1;metalNodes2;metalNodes3;metalNodes5]);
metalNodes = [];
%QMnodes = [semiNodes3];
%sQMnodes = [semiNodes10];
%sQMnodes1= [semiNodes12;semiNodes13];
semiintf   = intersect(semiNodesn,semiNodesp);

%semiNodes  = unique([semiNodesn;semiNodesp]);
%semiNodes  = [];

issemiintf = false(Nnode,1);
issemiintf(semiintf) = true;

doprf  = zeros(Nnode,1);

for iii=1:length(semiintf)
    ajlk = nodeLinks{semiintf(iii)}(1,:);
    ajnd = nodeLinks{semiintf(iii)}(2,:);
    for ttt=1:length(ajnd)
       if issemiintf(ajnd(ttt)) == false
          disif  = linkL(ajlk(ttt));
          doprf(ajnd(ttt)) = doprf(ajnd(ttt)) + exp(-(disif/0.2)^2);
       end
    end
end

QMnodes = [semiNodes3];
sQMnodes = [semiNodes10];

interfNodes = intersect(semiNodes,metalNodes); % interface between metal and semiconductor
dirSemiNodes = interfNodes;
isDirSemiNodes = false(Nnode,1);
isDirSemiNodes(dirSemiNodes) = true;

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
    if (isSqmnodes(n1)+isSqmnodes(n2)) > 0
      if (isSqmnodes(n1)+isSqmnodes(n2))==1
        sQMlinks = [sQMlinks;ilk];
      end
      QMlinks = [QMlinks;ilk];
      isSQMlinks(ilk)  = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idSemiNodes = false(Nnode,1);
idSemiNodes(semiNodes) = true;
idMetalNodes = false(Nnode,1);
idMetalNodes(metalNodes) = true;

%%%% define doping profile %%%%%%%%%%%%%%
%doping = dopeUniform(semiNodesn,semiNodesp,semiintf,doprf); % only uniform doping
%doping = dopeUniform(semiNodes,QMnodes);
doping = dopeUniformnp(contacts{2}.nodes,contacts{1}.nodes,semiintf);

% eqnNodes = setdiff((1:Nnode)',dirNodes);
% eqnSemiNodes = setdiff(semiNodes,dirSemiNodes);

bndLinks = [];
for i = 1:Nlink
    n1 = links(i,1); n2 = links(i,2);
    if isBndNodes(n1) && isBndNodes(n2)
       bndLinks = [bndLinks;i]; 
    end
end

%%% material types of each volume (semiconductor = 1, metal = 2, insulator =3 %%%
volumeM = 3*ones(Nvolume,1);
for i = 1:Nvolume
   vNodes = volumeNodes(i,:);
   if all(idSemiNodes(vNodes))
      volumeM(i) = 1; 
   end
   if all(idMetalNodes(vNodes))
      volumeM(i) = 2; 
   end
end

isqmvolm = false(Nvolume,1);
for i = 1:Nvolume
   vNodes = volumeNodes(i,:);
   if all(istqmnodes(vNodes))
      isqmvolm(i) = true;
   end
end

issemvolm_p = false(Nvolume,1);
issemvolm_n = false(Nvolume,1);


idseminodevolmn = false(Nvolume,1);
idseminodevolmn(seminodevolmn) = true;
idseminodevolmp = false(Nvolume,1);
idseminodevolmp(seminodevolmp) = true;

seminodevolm=[];

for i = 1:Nvolume
   vNodes = volumeNodes(i,:);
   xc=sum(nodes(vNodes,1))/8;
   yc=sum(nodes(vNodes,2))/8;
   zc=sum(nodes(vNodes,3))/8;

   if idseminodevolmp(i)
      seminodevolm = [seminodevolm;[xc,yc,zc,0]];
      issemvolm_p(i) = true;
   end
   if idseminodevolmn(i)
      seminodevolm = [seminodevolm;[xc,yc,zc,1]];
      issemvolm_n(i) = true;
   end
end


display('End solver initialization');
