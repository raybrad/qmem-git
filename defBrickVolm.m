function brikvolm = defBrickNodes(p1,p2)

%%% provide all the nodes within a cube defined by its lower left and upper right
%%% points

global kx ky kz;
global nodes links contacts;
global Nnode Nlink Nvolume;
global nodeM volumeM;
global bndNodes intNodes dirNodes eqnNodes edgeNodes;
global eqnLinks  bndLinks ;
global isBndNodes isDirNodes dcVolDirNodes acVolDirNodes nodeLinks volumeNodes;
global doping;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;

global nodeLinks;
global linkL;

x1 = p1(1); y1 = p1(2); z1 = p1(3);
x2 = p2(1); y2 = p2(2); z2 = p2(3);

brknodes = [];
for z = z1:z2
   for y = y1:y2
      for x = x1:x2
         brknodes = [brknodes;co2id([x,y,z])];
      end
   end
end

brikvolm=[];

for n1 = brknodes.'
   ajvol_n1 = nodeVolumes{n1}(1,:);

   for j=1:length(ajvol_n1)
     xc=0; yc=0; zc=0;

     vNodes = volumeNodes(ajvol_n1(j),:);

     for k=1:length(vNodes)
       [nx,ny,nz]=id2co(vNodes(k));
       xc=xc+nx;
       yc=yc+ny;
       zc=zc+nz;;
     end

     xc=xc/8;
     yc=yc/8;
     zc=zc/8;
   
     if xc>x1 && xc<x2 && yc>y1 && yc<y2 && zc>z1 && zc<z2
      brikvolm= unique([brikvolm;ajvol_n1(j)]);   
     end
   end
end


