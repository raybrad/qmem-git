function [Nodex,Nodey,Nodez] = RealToNodes(RealCoor,x_coor,y_coor,z_coor)

%%Tranform a real position into a node position	
%%RealCoor: the Real coordinate
realx = RealCoor(1); realy = RealCoor(2); realz = RealCoor(3);

nx=length(x_coor);
ny=length(y_coor);
nz=length(z_coor);
for x=1:nx-1
  if (realx>=x_coor(x) && realx< x_coor(x+1) )
	Nodex=x;
  elseif(realx==x_coor(x+1))
	Nodex=x+1;
  end
end  

for y=1:ny-1
  if (realy>=y_coor(y) && realy< y_coor(y+1) )
	Nodey=y;
  elseif(realy==y_coor(y+1))
	Nodey=y+1;
  end
end  

for z=1:nz-1
  if (realz>=z_coor(z) && realz< z_coor(z+1) )
	Nodez=z;
  elseif(realz==z_coor(z+1))
	Nodez=z+1;
  end
end 
