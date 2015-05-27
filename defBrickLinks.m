function inlinks = defBrickLinks(p1,p2)

%%% provide all the nodes within a cube defined by its lower left and upper right
%%% points

global nodes links;
global Nnode Nlink;

innodes   = defBrickNodes(p1,p2);
isinnodes         = false(Nnode,1);
isinnodes(innodes)= true; 

inlinks = [];


for i = 1:Nlink
    if isinnodes(links(i,1)) && isinnodes(links(i,2))
       inlinks=[inlinks;i];
    end
end
