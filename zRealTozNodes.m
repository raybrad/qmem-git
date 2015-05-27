function Nodez = zRealTozNodes(realz,z_coor)

%%Tranform a real position into a node position	
%%RealCoor: the Real coordinate

nz=length(z_coor);

for z=1:nz-1
  if (realz>=z_coor(z) && realz< z_coor(z+1) )
	Nodez=z;
  elseif(realz==z_coor(z+1))
	Nodez=z+1;
  end
end 
