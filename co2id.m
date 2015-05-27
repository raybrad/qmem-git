function nid = co2id(nco)

global kx ky kz;
nx = nco(:,1); ny = nco(:,2); nz = nco(:,3);
nid = nx+(ny-1)*(kx+1)+(nz-1)*(kx+1)*(ky+1);

