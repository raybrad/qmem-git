function corpedGeometry(periodx, periody, periodz)

display('Start geometry correction for period condition');

global kx ky kz;
global nodes links;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks volumeSurfs linkVolumes nodeVolumes;
global Nnode Nlink Nsurf Nvolume;
global nodeV linkL linkS dlinkL linkCenter surfCenter;

global bndNodes edgeNodes dirNodes;
global isBndNodes isDirNodes dcVolDirNodes acVolDirNodes;
global doping;
global semiNodes;
global dirSemiNodes isDirSemiNodes;
global bndLinks;
global metalNodes;


numx=kx+1;
numy=ky+1;
numz=kz+1;

bboundary   =[];
lboundary   =[];
zboundary   =[];
pedbndnodes =[];

if periodz
   bottom_boundary = 1:numx*numy;
   top_boundary    = (numx*numy*(numz-1))+1:numx*numy*numz;
   zboundary       = [bottom_boundary top_boundary ];
   zboundary       = unique(zboundary);
   zboundary       = zboundary';
end  

if periodx
   left_boundary   = 1:numx:(numx*numy*numz-numx+1);
   right_boundary  = numx:numx:numx*numy*numz;
   xboundary       = [left_boundary right_boundary ];
   xboundary       = unique(xboundary);
   xboundary       = xboundary';
end

if periody
 for k = 1:numy
    for j = 1:numz
        for i = 1:numx
            dom(i,j,k) = i+ numx*(k-1) + (j-1)*numx*numy;
        end
    end
 end
 front_boundary    = reshape(dom(:,:,1),numx*numz,1)';
 back_boundary     = reshape(dom(:,:,numy),numx*numz,1)';
   yboundary       = [front_boundary back_boundary ];
   yboundary       = unique(yboundary);
   yboundary       = yboundary';
end

if periodx + periody + periodz ==2
   xycorss    = intersect(xboundary,yboundary);
   yzcorss    = intersect(yboundary,zboundary);
   zxcorss    = intersect(zboundary,xboundary);
   pedbndnodes= unique([xycorss;yzcorss;zxcorss]);
end

xboundaryn = setdiff(xboundary,pedbndnodes);
yboundaryn = setdiff(yboundary,pedbndnodes);
zboundaryn = setdiff(zboundary,pedbndnodes);

for n1 = xboundaryn.'
       ajlk_n1 = nodeLinks{n1}(1,:);
       ajnd_n1 = nodeLinks{n1}(2,:);
       ajvol_n1 = nodeVolumes{n1}(1,:);
       ajvolV_n1 = nodeVolumes{n1}(2,:);
end

for n1 = yboundaryn.'
        
       [nx,ny,nz]=id2co(n1);
       
       if ny == 1
          nl1 = nodeLinks{n1};
          nv1 = nodeVolumes{n1};
          n2=co2id([nx,numy,nz]);
          nl2 = nodeLinks{n2};
          nv2 = nodeVolumes{n2};
         
          nlt = [nl1,nl2];
          nvt = [nv1,nv2];
             
          nodeLinks{n1}   = nlt;
          nodeLinks{n2}   = nlt;          
          nodeVolumes{n1} = nvt;
          nodeVolumes{n2} = nvt;
          

          linkSurfs
          
          linkVolumes


       end
end


display('End geometry correction for period condition');

