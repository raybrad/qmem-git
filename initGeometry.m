function initGeometry()

display('Start geometry initialization');

global kx ky kz;
global nodes links;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks volumeSurfs linkVolS nodeVolV;
global Nnode Nlink Nsurf Nvolume;
global nodeV linkL linkS dlinkL linkCenter surfCenter;
%%%new global variables for light source and bnd
global dirxLinks diryLinks dirzLinks;
global EsurfLinks BsurfLinks;
global lightdirection;

tStart=tic;
tic;
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
case 'kzEyBx'	
Leftsurf= defBrickLinks([1,1,1],[kx+1,1,kz+1]);
Rightsurf= defBrickLinks([1,ky+1,1],[kx+1,ky+1,kz+1]);
EsurfLinks=[Leftsurf;Rightsurf];	%E

Backsurf=defBrickLinks([1,1,1],[1,ky+1,kz+1]);
Frontsurf=defBrickLinks([kx+1,1,1],[kx+1,ky+1,kz+1]);
BsurfLinks=[Backsurf;Frontsurf];	%B
case 'kzExBy'	
Leftsurf= defBrickLinks([1,1,1],[kx+1,1,kz+1]);
Rightsurf= defBrickLinks([1,ky+1,1],[kx+1,ky+1,kz+1]);
BsurfLinks=[Leftsurf;Rightsurf];	%B

Backsurf=defBrickLinks([1,1,1],[1,ky+1,kz+1]);
Frontsurf=defBrickLinks([kx+1,1,1],[kx+1,ky+1,kz+1]);
EsurfLinks=[Backsurf;Frontsurf];	%E
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adjacent links and nodes of nodes               connectlinks   connectnodes                    
nodeLinks = cell(Nnode,1); %format:{target node} [link1,link2...;node1,node2...],how each node connects to other nodes

% centers of links
linkCenter = zeros(Nlink,3);
for i = 1:Nlink
    n1 = links(i,1); n2 = links(i,2);
    nodeLinks{n1} = [nodeLinks{n1},[i;n2]];
    nodeLinks{n2} = [nodeLinks{n2},[i;n1]];
    linkCenter(i,:) = (nodes(n1,:)+nodes(n2,:))/2;
    links(i,3) = find(nodes(n1,:)-nodes(n2,:));
end

display('End nodeLinks');
%node connectivity of surfaces
surfNodes = zeros(Nsurf,4);
surfCenter = zeros(Nsurf,3);
for i = 1:Nsurf
    i1 = (i-1)*4+1;
    lk = surfLinks(i1:i1+3,2);
    nd = links(lk,1:2);
    surfNodes(i,1:2) = nd(1,:);
    surfNodes(i,3) = setdiff(nd(2,:),nd(1,:));
    surfNodes(i,4) = setdiff(nd(3,:),surfNodes(i,3));
    surfCenter(i,:) = sum(nodes(surfNodes(i,:),:))/4;
end
display('End surfNodes');

%nodes(Nnodes,xyz coord)
%length of links	 (links(Nlink,sta),
linkL = sqrt(sum((nodes(links(:,1),:)-nodes(links(:,2),:)).^2,2));       %linkL(Nlink)

%link connectivity of surfaces
surfLinks = reshape(surfLinks(:,2),4,Nsurf)';

display(['time for part 1:']);
toc;
tic;

%link connectivity of volumes & dual area of links
nodeV = zeros(Nnode,1);
linkS = zeros(Nlink,1);
dlinkL = zeros(Nsurf,1);
%nodeVolV = zeros(Nnode,Nvolume);
ntripletsNV=8*Nvolume;
rowNV=zeros(ntripletsNV,1);
colNV=zeros(ntripletsNV,1);
valNV=zeros(ntripletsNV,1);
ntripletsNV=0;

for i = 1:size(volumeNodes,1)
   nodeV(volumeNodes(i,2)) = nodeV(volumeNodes(i,2))+volumeNodes(i,3);	% 1/8+1/8+... 
%   nodeVolV(volumeNodes(i,2),volumeNodes(i,1)) = volumeNodes(i,3);    %1/8
   ntripletsNV=ntripletsNV+1;
   rowNV(ntripletsNV)=volumeNodes(i,2);
   colNV(ntripletsNV)=volumeNodes(i,1);
   valNV(ntripletsNV)=volumeNodes(i,3);
end
nodeVolV=sparse(rowNV(1:ntripletsNV),colNV(1:ntripletsNV),valNV(1:ntripletsNV),Nnode,Nvolume);
display('finish generating nodeVolV');
display(['time for part 2:']);
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
%  linkVolS(volumeLinks(i,2),volumeLinks(i,1)) = volumeLinks(i,3);
   	    %(link,volumeid)                      surface 1/4
   ntripletslVS=ntripletslVS+1;
   rowlVS(ntripletslVS)=volumeLinks(i,2);	%link
   collVS(ntripletslVS)=volumeLinks(i,1);	%
   vallVS(ntripletslVS)=volumeLinks(i,3);
end
linkVolS=sparse(rowlVS(1:ntripletslVS),collVS(1:ntripletslVS),vallVS(1:ntripletslVS),Nlink,Nvolume);
display('finish generating linkVolS');

display(['time for part 3:']);
toc;
tic;

for i = 1:size(volumeSurfs,1)
%which cube%volumeSurf~ one half of the prependicular length  in the corresponding cube with respect to a surf
   dlinkL(volumeSurfs(i,2)) = dlinkL(volumeSurfs(i,2))+volumeSurfs(i,3); % like the average distance of neighbor cube  
end
tmp = volumeNodes(:,2);
volumeNodes = reshape(tmp,8,Nvolume)';
tmp = volumeLinks(:,2);
volumeLinks = reshape(tmp,12,Nvolume)';

% surface connectivity of volumes & length of dual link through the centers
% of surfaces
tmp = volumeSurfs(:,2);
volumeSurfs = reshape(tmp,6,Nvolume)';
											           %the other link on that surf
% adjacent surfaces of links. Format: [surf;otherlinks]   ~  linkSurfs{link} ~{surf index the link belong to;link1;link2;link3}
							%		      {surf index the link belong to;link1;link2;link3}	
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
display(['time for part 4:']);
toc;
%tic;
%%% adjacent volumes of nodes and links. Format: [volumeid;associated volume or area]
%nodeVolumes = cell(Nnode,1);
%linkVolumes = cell(Nlink,1);
%for i=1:Nnode
%IndexNodes=find(nodeVolV(i,:));
%nodeVolmes{i}=zeros(2,length(IndexNodes));
%nodeVolumes{i}(1,:)=IndexNodes;
%nodeVolumes{i}(2,:)=nodeVolV(i,IndexNodes);
%end

%for i=1:Nlink
%IndexLinks=find(linkVolS(i,:));
%linkVolmes{i}=zeros(2,length(IndexLinks));
%linkVolumes{i}(1,:)=IndexLinks;
%linkVolumes{i}(2,:)=linkVolS(i,IndexLinks);
%end


%linkVolumes = cell(Nlink,1);
%nodeVolumes = cell(Nnode,1);
%for i = 1:Nvolume
%    nd = volumeNodes(i,:);
%   for j = 1:length(nd)
%       nodeVolumes{nd(j)} = [nodeVolumes{nd(j)},[i;nodeVolV(nd(j),i)]]; % volumeid1   1   2   ''''''' 
%   end                                                                  % volume1     V/8 V/8 '''''''''
%   lk = volumeLinks(i,:);                                                      
%  for j = 1:length(lk)
%      linkVolumes{lk(j)} = [linkVolumes{lk(j)},[i;linkVolS(lk(j),i)]];
%      			%	{lk}	% volumeid1    (1,1)	volumeid2     (1,2)    (1,3)...
%      		 	% usual 2*4	% surface area1(2,1)	surface area2 (2,2)    (2,3)...
%  end
%end

%display('End nodeVolumes & linkVolumes');
%display(['time for part 5:']);
%toc;
%tic;


%%% scale down to nm range %%%%%%

display('End geometry initialization');
toc(tStart);
