function  dembv=tdembv(VT)

global kx ky kz;
global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;

nst=length(VT(1,:));
dembv=[];

for i=1:nst
  V=VT(:,i);
  input_vs=reshape(V,kx+1,ky+1,kz+1);
  qm_vs = input_vs(qmx1:qmx2,qmy1:qmy2,qmz1:qmz2);
  leftnodes=defBrickNodes([qmx1,qmy1,qmz1],[qmx1,qmy2,qmz2]);
  rghtnodes=defBrickNodes([qmx2,qmy1,qmz1],[qmx2,qmy2,qmz2]);
  ltv=mean(V(leftnodes));
  rtv=mean(V(rghtnodes));
  dtv=ltv-rtv;
  dembv=[dembv;dtv];
end
