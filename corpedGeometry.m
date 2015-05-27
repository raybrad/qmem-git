function corpedGeometry(nperiodd)

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

ybndnodes   = defBrickNodes([2,1,2],[kx,1,kz]);
zbndnodes   = defBrickNodes([2,2,1],[kx,ky,1]);
yzbndnodes  = defBrickNodes([2,1,1],[kx,1,1]);

lnboundary  =[];

zboundarynl =[];
xboundarynl =[];
yboundarynl =[];


if periodz
   bboundary = 1:numx*numy;
end  

if periodx
   lboundary   = 1:numx:(numx*numy*numz-numx+1);
end

if periody
 for k = 1:numy
    for j = 1:numz
        for i = 1:numx
            dom(i,j,k) = i+ numx*(k-1) + (j-1)*numx*numy;
        end
    end
 end
 fboundary    = reshape(dom(:,:,1),numx*numz,1)';
end

if periodx + periody + periodz ==2
   xycorss    = intersect(lboundary,fboundary);
   yzcorss    = intersect(bboundary,fboundary);
   zxcorss    = intersect(bboundary,lboundary);
   lnboundary = unique([xycorss;yzcorss;zxcorss]);
end

xboundarynl = setdiff(lboundary,edgeNodes);
yboundarynl = setdiff(fboundary,edgeNodes);
zboundarynl = setdiff(bboundary,edgeNodes);

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

