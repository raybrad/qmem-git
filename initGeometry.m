function initGeometry()

display('--------------------------------------------------------------------');
display('Start geometry initialization');
tStart=tic;
tic;

global kx ky kz;
global nodes links;
global nodeLinks linkSurfs surfLinks volumeNodes volumeLinks volumeSurfs linkVolumes nodeVolumes;
global Nnode Nlink Nsurf Nvolume;
global nodeV linkL linkS dlinkL;
%%%new global variables for light source and bnd
global dirxLinks diryLinks dirzLinks;
global EsurfLinks BsurfLinks;
global lightdirection;

nodes = sortrows(nodes);
nodes = nodes(:,2:4);
links = links(:,2:3);
Nnode = size(nodes,1);
Nlink = size(links,1);
Nsurf = size(surfLinks,1)/4;
Nvolume = size(volumeNodes,1)/8;
kx = length(unique(nodes(:,1)))-1;
ky = length(unique(nodes(:,2)))-1;
kz = length(unique(nodes(:,3)))-1;

%%%define directional Links,light sheet
dirxLinks=zeros(Nlink/3,1);
diryLinks=zeros(Nlink/3,1);
dirzLinks=zeros(Nlink/3,1);
    tmpcounter = 0;    
    indcounterx= 0;
    indcountery= 0;
    indcounterz= 0;
    for k = 1:kz+1
        for j = 1:ky+1
            for i = 1:kx+1
                if (i <=kx) 
                tmpcounter = 1 + tmpcounter;
		indcounterx=indcounterx+1;
                dirxLinks(indcounterx)=tmpcounter;
                end
                if (j <=ky) 
                tmpcounter = 1 + tmpcounter;
		indcountery=indcountery+1;
                diryLinks(indcountery)=tmpcounter;
                end
                if (k <=kz) 
                tmpcounter = 1 + tmpcounter;
		indcounterz=indcounterz+1;
                dirzLinks(indcounterz)=tmpcounter;
                end
            end
        end
   end
display('End directional links collection');
toc;
tic;
%%%%%%%%define E surface boundary and B surface boundary
switch lightdirection
case 'kxEyBz'
Leftsurf= defBrickLinks([1,1,1],[kx+1,1,kz+1]);
Rightsurf= defBrickLinks([1,ky+1,1],[kx+1,ky+1,kz+1]);
EsurfLinks=[Leftsurf;Rightsurf];	%E

Botsurf= defBrickLinks([1,1,1],[kx+1,ky+1,1]);
Topsurf= defBrickLinks([1,1,kz+1],[kx+1,ky+1,kz+1]);
BsurfLinks=[Botsurf;Topsurf];  %B
case 'kxEzBy'
Leftsurf= defBrickLinks([1,1,1],[kx+1,1,kz+1]);
Rightsurf= defBrickLinks([1,ky+1,1],[kx+1,ky+1,kz+1]);
BsurfLinks=[Leftsurf;Rightsurf];	%B

Botsurf= defBrickLinks([1,1,1],[kx+1,ky+1,1]);
Topsurf= defBrickLinks([1,1,kz+1],[kx+1,ky+1,kz+1]);
EsurfLinks=[Botsurf;Topsurf];  %E
case 'kyExBz'
Backsurf=defBrickLinks([1,1,1],[1,ky+1,kz+1]);
Frontsurf=defBrickLinks([kx+1,1,1],[kx+1,ky+1,kz+1]);
EsurfLinks=[Backsurf;Frontsurf];	%E

Botsurf= defBrickLinks([1,1,1],[kx+1,ky+1,1]);
Topsurf= defBrickLinks([1,1,kz+1],[kx+1,ky+1,kz+1]);
BsurfLinks=[Botsurf;Topsurf];  %B
case 'kyEzBx'
Backsurf=defBrickLinks([1,1,1],[1,ky+1,kz+1]);
Frontsurf=defBrickLinks([kx+1,1,1],[kx+1,ky+1,kz+1]);
BsurfLinks=[Backsurf;Frontsurf];	%B

Botsurf= defBrickLinks([1,1,1],[kx+1,ky+1,1]);
Topsurf= defBrickLinks([1,1,kz+1],[kx+1,ky+1,kz+1]);
EsurfLinks=[Botsurf;Topsurf];  %E
case {'kzEyBx','specialkzEyBx'}
Leftsurf= defBrickLinks([1,1,1],[kx+1,1,kz+1]);
Rightsurf= defBrickLinks([1,ky+1,1],[kx+1,ky+1,kz+1]);
EsurfLinks=[Leftsurf;Rightsurf];	%E

Backsurf=defBrickLinks([1,1,1],[1,ky+1,kz+1]);
Frontsurf=defBrickLinks([kx+1,1,1],[kx+1,ky+1,kz+1]);
BsurfLinks=[Backsurf;Frontsurf];	%B
case {'kzExBy','specialkzExBy'}	
Leftsurf= defBrickLinks([1,1,1],[kx+1,1,kz+1]);
Rightsurf= defBrickLinks([1,ky+1,1],[kx+1,ky+1,kz+1]);
BsurfLinks=[Leftsurf;Rightsurf];	%B

Backsurf=defBrickLinks([1,1,1],[1,ky+1,kz+1]);
Frontsurf=defBrickLinks([kx+1,1,1],[kx+1,ky+1,kz+1]);
EsurfLinks=[Backsurf;Frontsurf];	%E
case 'nolight'
    BsurfLinks=[];
    EsurfLinks=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('light source links');
toc;
tic;
% adjacent links and nodes of nodes               connectlinks   connectnodes                    
nodeLinks = cell(Nnode,1); %format:{target node} [link1,link2...;node1,node2...],how each node connects to other nodes
for i = 1:Nlink
    n1 = links(i,1); n2 = links(i,2);
    nodeLinks{n1} = [nodeLinks{n1},[i;n2]];
    nodeLinks{n2} = [nodeLinks{n2},[i;n1]];
    links(i,3) = find(nodes(n1,:)-nodes(n2,:));
end

display('End nodeLinks');
display(['time for nodeLinks:']);
toc;
tic;

%length of links	 (links(Nlink,sta),
linkL = sqrt(sum((nodes(links(:,1),:)-nodes(links(:,2),:)).^2,2));       %linkL(Nlink)

%link connectivity of surfaces
surfLinks = reshape(surfLinks(:,2),4,Nsurf)';

display(['time for surfLinks:']);
toc;
tic;

%link connectivity of volumes & dual area of links
nodeV = zeros(Nnode,1);
linkS = zeros(Nlink,1);
dlinkL = zeros(Nsurf,1);
ntripletsNV=8*Nvolume;
rowNV=zeros(ntripletsNV,1);
colNV=zeros(ntripletsNV,1);
valNV=zeros(ntripletsNV,1);
ntripletsNV=0;

for i = 1:size(volumeNodes,1)
   nodeV(volumeNodes(i,2)) = nodeV(volumeNodes(i,2))+volumeNodes(i,3);	% 1/8+1/8+... 
   ntripletsNV=ntripletsNV+1;
   rowNV(ntripletsNV)=volumeNodes(i,2);
   colNV(ntripletsNV)=volumeNodes(i,1);
   valNV(ntripletsNV)=volumeNodes(i,3);
end
nodeVolV=sparse2(rowNV(1:ntripletsNV),colNV(1:ntripletsNV),valNV(1:ntripletsNV),Nnode,Nvolume);
display(['time for nodeVolV:']);
toc;
tic;

%linkVolS = zeros(Nlink,Nvolume);
ntripletslVS=size(volumeLinks,1);
rowlVS=zeros(ntripletslVS,1);
collVS=zeros(ntripletslVS,1);
vallVS=zeros(ntripletslVS,1);
ntripletslVS=0;
for i = 1:size(volumeLinks,1)
%one fourth of the area of the orthogonal surface in the corresponding cube with respect to link
   linkS(volumeLinks(i,2)) = linkS(volumeLinks(i,2))+volumeLinks(i,3); 
   ntripletslVS=ntripletslVS+1;
   rowlVS(ntripletslVS)=volumeLinks(i,2);	%link
   collVS(ntripletslVS)=volumeLinks(i,1);	%volumeid
   vallVS(ntripletslVS)=volumeLinks(i,3);
end
linkVolS=sparse2(rowlVS(1:ntripletslVS),collVS(1:ntripletslVS),vallVS(1:ntripletslVS),Nlink,Nvolume);
display(['time for linkVolS:']);
toc;
tic;

%which cube%volumeSurf~ one half of the prependicular length  in the corresponding cube with respect to a surf
for i = 1:size(volumeSurfs,1)
   dlinkL(volumeSurfs(i,2)) = dlinkL(volumeSurfs(i,2))+volumeSurfs(i,3); % like the average distance of neighbor cube  
end
tmp = volumeNodes(:,2);
volumeNodes = reshape(tmp,8,Nvolume)';
tmp = volumeLinks(:,2);
volumeLinks = reshape(tmp,12,Nvolume)';
tmp = volumeSurfs(:,2);
volumeSurfs = reshape(tmp,6,Nvolume)';
display(['time for dLinkL:']);
toc;
tic;
									      % the other link on that surf
% adjacent surfaces of links. Format: [surf;otherlinks]   ~  linkSurfs{link} ~{surf index the link belong to;link1;link2;link3}
								              %{surf index the link belong to;link1;link2;link3}	
linkSurfs = cell(Nlink,1);
for i = 1:Nsurf
    l1 = surfLinks(i,1); l2 = surfLinks(i,2); l3 = surfLinks(i,3); l4 = surfLinks(i,4);
    if any(links(l1,2) == links(l2,:))
        linkSurfs{l1} = [linkSurfs{l1},[i;l2;l3;l4]];
    else
        linkSurfs{l1} = [linkSurfs{l1},[i;l4;l3;l2]];
    end
    if any(links(l2,2) == links(l3,:))
        linkSurfs{l2} = [linkSurfs{l2},[i;l3;l4;l1]];
    else
        linkSurfs{l2} = [linkSurfs{l2},[i;l1;l4;l3]];
    end
    if any(links(l3,2) == links(l4,:))
        linkSurfs{l3} = [linkSurfs{l3},[i;l4;l1;l2]];
    else
        linkSurfs{l3} = [linkSurfs{l3},[i;l2;l1;l4]];
    end
    if any(links(l4,2) == links(l1,:))
        linkSurfs{l4} = [linkSurfs{l4},[i;l1;l2;l3]];
    else
        linkSurfs{l4} = [linkSurfs{l4},[i;l3;l2;l1]];
    end
end
display('End linkSurfs');
toc;
tic;
%%%%%
linkVolumes = cell(Nlink,1);
nodeVolumes = cell(Nnode,1);
for i = 1:Nvolume*8
        nodeVolumes{rowNV(i)} = [nodeVolumes{rowNV(i)},[colNV(i);valNV(i)]]; % 1   1   1  1 1 1 1 1
end                                                                             % V/8 V/8 '''''''''
for i=1:Nvolume*12
        linkVolumes{rowlVS(i)}= [linkVolumes{rowlVS(i)},[collVS(i);vallVS(i)]];
						% volumeid;surface area
end
display('End linkVolumes nodeVolmes');
toc;
tic;
%%% scale down to nm range %%%%%%

display('End geometry initialization');
toc(tStart);
display('--------------------------------------------------------------------');
