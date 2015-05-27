function [nx,ny,nz] = id2co(nid)

global kx ky kz;
global Node;
if nid==Node
   nx = kx+1;
   ny = ky+1;
   nz = kz+1;
else
   j0 = floor((nid-1)/((kx+1)*(ky+1)));
   nz = j0+1;
   k0 = nid - j0*(kx+1)*(ky+1);
   l0 = floor((k0-1)/(kx+1));
   ny = l0+1;
   nx = k0-l0*(kx+1);
end
