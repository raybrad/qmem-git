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
global metalNodes;

%%%
global dirxLinks diryLinks dirzLinks;
global JLinks J_amp;
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
%R
X1 = 6; X2 = 11; Y1 = 6; Y2 = 11; Z1 = 6; Z2 = 11;Xms=3;Xme=14;
qmx1 = X1; qmx2 = X2; qmy1 = Y1; qmy2 = Y2; qmz1 = Z1; qmz2 = Z2;
Radius=2.0;
Height=4.0;
Xcen=8;Ycen=8;Zcen=8;
avoltage   = 0.0;
avoltageac = 0.0;
%R: set the boundary voltage to 0 now, for light
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
veclayer= defBrickLinks([1,1,kz+1],[kx+1,ky+1,kz+1]);
JLinks=intersect(veclayer,diryLinks);
J_amp=1e6*ones(length(JLinks),1)/scl.s_J;
%%%%%%%%%%%

%metalNodes = defBrickNodes([X1,Y1,Z1],[X2,Y2,Z2]); % metal from contact to semi surface
metalNodes=defCylinderNodes([Xcen,Ycen,Zcen],Radius,Height,'z');

semiNodes  = [];
semiintf =[];

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

QMnodes = [];
sQMnodes = [];

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
idSemiNodes = false(Nnode,1);
idSemiNodes(semiNodes) = true;
idMetalNodes = false(Nnode,1);
idMetalNodes(metalNodes) = true;

%%%% define doping profile %%%%%%%%%%%%%%
%doping = dopeUniform(semiNodesn,semiNodesp,semiintf,doprf); % only uniform doping
doping = dopeUniform(semiNodes,QMnodes);

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

display('End solver initialization');
