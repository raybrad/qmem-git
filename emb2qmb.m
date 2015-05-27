function  emb2qmb(V,effcd,currtime,dt)

global kx ky kz;
global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;
%
global qkx qky qkz;
global ingvl ingvr ingvt;

ix = qmx2-qmx1+1; iy = qmy2-qmy1+1; iz = qmz2-qmz1+1;


%qmx = 4; qmy = 1; qmz=1;
qmx = qmx2-qmx1; qmy = qmy2-qmy1; qmz = qmz2-qmz1;

  input_vs=reshape(V,kx+1,ky+1,kz+1);
  qm_vs = input_vs(qmx1:qmx2,qmy1:qmy2,qmz1:qmz2);

  input_cd=reshape(effcd,kx+1,ky+1,kz+1);
  qm_cd = input_cd(qmx1:qmx2,qmy1:qmy2,qmz1:qmz2);

  leftnodes=defBrickNodes([qmx1,qmy1,qmz1],[qmx1,qmy2,qmz2]);
  rghtnodes=defBrickNodes([qmx2,qmy1,qmz1],[qmx2,qmy2,qmz2]);
  topnodes =defBrickNodes([qmx1,qmy1,qmz2],[qmx2,qmy2,qmz2]);

  ltv=mean(V(leftnodes));

  rtv=mean(V(rghtnodes));

  ttv=mean(V(topnodes));

  if currtime > 0
   ingvl = ingvl + ltv * dt;
   ingvr = ingvr + rtv * dt;
   ingvt = ingvt + ttv * dt;
  end

  [xi,yi,zi] = meshgrid(0:qmy/(qky-1):qmy,0:qmx/(qkx-1):qmx,0:qmz/(qkz-1):qmz);
  [x,y,z]=meshgrid(0:qmy/(iy-1):qmy,0:qmx/(ix-1):qmx,0:qmz/(iz-1):qmz);

  emboundqm3d=interp3(x,y,z,qm_vs,xi,yi,zi);
  embound=reshape(emboundqm3d,qkx*qky*qkz,1);

  emboundqm3dcd=interp3(x,y,z,qm_cd,xi,yi,zi);
  emboundcd=reshape(emboundqm3dcd,qkx*qky*qkz,1);

  embound=[embound,emboundcd];
 
%  [xi,yi,zi] = meshgrid(0:qmx/qkx:qmx, 0:qmy/qky:qmy, 0:qmz/qkz:qmz);

%  [x,y,z]=meshgrid(0:0.5:qmx,0:0.25:qmy,0:0.25:qmz);

%  embound1=interp3(x,y,z,qm_vs,xi,yi,zi);
%  embound=reshape(qm_vs,ix*iy*iz,1);
  fp = fopen('embound-v','w');
  fprintf(fp,'%f\n',currtime);
  fprintf(fp,'%f\n',dt);
  fprintf(fp,'%13.5e %13.5e %13.5e \n',ltv,rtv,ttv);
  fprintf(fp,'%13.5e %13.5e %13.5e \n',ingvl,ingvr,ingvt);
  fclose(fp);

  fp = fopen('embound1','w');
  fprintf(fp,'%d %d %d\n',qkx,qky,qkz);
  fclose(fp);
  save('embound2','embound','-ASCII','-double');
  ! cat embound1 embound2 > embound;
  ! rm embound1 embound2;

  fp=fopen('voltage.dat','a+');
  fprintf(fp,'%f %13.5e %13.5e %13.5e \n',currtime,ltv,rtv,ttv);
  fclose(fp);

