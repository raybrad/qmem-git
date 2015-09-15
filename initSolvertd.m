function initSolvertd()
display('--------------------------------------------------------------------');
display('Start solver initialization');
tStart=tic;
tic;
%%%%%%%%%%%%%%%%%%Global var%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global kx ky kz scl;
global dt;
global nodes links contacts;
global Nnode Nlink Nvolume; 
global volumeM;
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
global qmRegionX1 qmRegionX2 qmRegionY1 qmRegionY2 qmRegionZ1 qmRegionZ2;
global isqmvolm;
global metalNodes metalLinks;

global x_coor y_coor z_coor;
global avoltage avoltageac;
global dirxLinks diryLinks dirzLinks;
global JLinks J_amp amplitude;
global lightdirection Posz;
global taskOpt;
global omega_p gamma_p;
global lepsr_1 lomega_1 lgamma_1 lepsr_2 lomega_2 lgamma_2;
global CJ01 CJ02 CJ11 CJ12 CJ13 CJ21 CJ22 CJ23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% define boundary and edge nodes %%%%%%%%%%%%%%%%%%%%%
[bndNodes,intNodes] = defBoundaryNodes(kx+1,ky+1,kz+1);	%bndNodes on 6 faces of the whole cube
isBndNodes = false(Nnode,1);
isBndNodes(bndNodes) = true;
edgeNodes = zeros(Nnode,1);
tmpcounter=0;
for i = bndNodes'
        ajnd = nodeLinks{i}(2,:); %connecting nodes to bndNodes
    if all(isBndNodes(ajnd))
    	tmpcounter=tmpcounter+1;
        edgeNodes(tmpcounter) = i; %edge nodes on the 12 edge lines
    end
end
edgeNodes(tmpcounter+1:Nnode)=[];
display('time for bndNodes edgeNodes');
toc;
tic;

%%%%%%%%%% structure definition (in terms of node index) %%%%%%%%%%%%%%%%%%

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
display('time for structrure and ac dc Nodes');
toc;
tic;
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
    fprintf('Position of light %d \n',tmpz);
	veclayer= defBrickLinks([1,1,tmpz],[kx+1,ky+1,tmpz]);
	JLinks=intersect(veclayer,dirxLinks);
case 'specialkzEyBx'	
	tmpz=zRealTozNodes(Posz,z_coor);
    fprintf('Position of light %d \n',tmpz);
	veclayer= defBrickLinks([1,1,tmpz],[kx+1,ky+1,tmpz]);
	JLinks=intersect(veclayer,diryLinks);
end
J_amp=ones(length(JLinks),1)/scl.s_J*amplitude;


%constant coefficient for plasmonic metal
CJ01=(2-gamma_p*dt)/(2+gamma_p*dt);
CJ02=dt*omega_p*omega_p/(2+gamma_p*dt);

CJ11=(2-lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt);
CJ12=(lgamma_1*dt-1)/(1+lgamma_1*dt);
CJ13=(lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt);

CJ21=(2-lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt);
CJ22=(lgamma_2*dt-1)/(1+lgamma_2*dt);
CJ23=(lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt);
%%%%%

if(taskOpt==1)
QMnodes = [];
sQMnodes = [];
else
QMnodes =defBrickNodes([qmx1,qmy1,qmz1],[qmx2,qmy2,qmz2]);
sQMnodes =defBrickNodes([qmx1+1,qmy1+1,qmz1+1],[qmx2-1,qmy2-1,qmz2-1]);
end

istqmnodes            = false(Nnode,1);
istqmnodes(sQMnodes)   = true;

isSqmnodes            = false(Nnode,1);
isSqmnodes(sQMnodes)  = true;

isqmnodes             = false(Nnode,1);
isqmnodes(sQMnodes)   = true;
isSQMlinks            = zeros(Nlink,1);

QMlinks   =[];
sQMlinks  =[];
currdlink = zeros(Nlink,6);

display('time for qm Nodes');
toc;
tic;

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
        if nx < ( qmRegionX1 + qmRegionX2 ) /2
         isSQMlinks(ilk)  = 1;
        else
         isSQMlinks(ilk)  = 2;
        end
     elseif n3==2 
        if ny < ( qmRegionY1 + qmRegionY2 )/2
         isSQMlinks(ilk)  = 3;        
        else
         isSQMlinks(ilk)  = 4;
        end
     elseif n3==3 
        if nz < ( qmRegionZ1 + qmRegionZ2 ) /2
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

display('time for qm Links');
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idMetalNodes = false(Nnode,1);
idMetalNodes(metalNodes) = true;
%idInsNodes = false(Nnode,1);
%idInsNodes(insNodes) = true;

%%%% define doping profile %%%%%%%%%%%%%%

%%%%%%%%%%%%
bndLinks = zeros(Nlink,1);
tmpcounter=0;
for i = 1:Nlink
    n1 = links(i,1); n2 = links(i,2);
    if isBndNodes(n1) && isBndNodes(n2)
       tmpcounter=tmpcounter+1;
       bndLinks(tmpcounter) = i; 
    end
end
bndLinks(tmpcounter+1:Nlink)=[];

metalLinks=zeros(Nlink,1);
tmpcounter=0;
for i = 1:Nlink
    n1 = links(i,1); n2 = links(i,2);
    if idMetalNodes(n1) && idMetalNodes(n2)
       tmpcounter=tmpcounter+1;
       metalLinks(tmpcounter) = i; 
    end
end
metalLinks(tmpcounter+1:Nlink)=[];

display('time for bndLinks metalLinks');
toc;
tic;
%%% material types of each volume %%%
volumeM = 2*ones(Nvolume,1);	%insulator 2 metal 1
for i = 1:Nvolume
   vNodes = volumeNodes(i,:);
   if all(idMetalNodes(vNodes))
      volumeM(i) = 1; 
   end
   if (taskOpt==2 && all(istqmnodes(vNodes)))
      volumeM(i) = 3; 
   end
end

isqmvolm = false(Nvolume,1);
for i = 1:Nvolume
   vNodes = volumeNodes(i,:);
   if all(istqmnodes(vNodes))
      isqmvolm(i) = true;
   end
end
display('time for volumeM');
toc;

display('End solver initialization');
toc(tStart);
display('--------------------------------------------------------------------');
