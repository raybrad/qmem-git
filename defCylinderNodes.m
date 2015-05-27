function nodesc = defCylinderNodes(p1,dd,hh,direction)

%%define the nodes of a cylinder with axis in z direction
%%p1: the center coordinate of cylinder, index
%%dd: the radius 
%%hh: the height

global kx ky kz;
global nodes;
x1 = p1(1); y1 = p1(2); z1 = p1(3);
tmpcounter=(kx+1)*(ky+1)*(z1-1)+(kx+1)*(y1-1)+x1;
x1_coor=nodes(tmpcounter,1);
y1_coor=nodes(tmpcounter,2);
z1_coor=nodes(tmpcounter,3);


nodesc = [];
tmpcounter=0;

switch direction

case 'x'
for z = 1:kz+1
   for y = 1:ky+1
      for x = 1:kx+1
         tmpcounter = 1 + tmpcounter;
	 x0=nodes(tmpcounter,1);
	 y0=nodes(tmpcounter,2);
	 z0=nodes(tmpcounter,3);
         ds=sqrt((z0-z1_coor)^2+(y0-y1_coor)^2);
         if (ds<=dd)&&(x0<=x1_coor+hh/2)&&(x0>=x1_coor-hh/2)
          nodesc = [nodesc;co2id([x,y,z])]; 
         end
      end
   end
end

case 'y'
for z = 1:kz+1
   for y = 1:ky+1
      for x = 1:kx+1
         tmpcounter = 1 + tmpcounter;
	 x0=nodes(tmpcounter,1);
	 y0=nodes(tmpcounter,2);
	 z0=nodes(tmpcounter,3);
         ds=sqrt((z0-z1_coor)^2+(x0-x1_coor)^2);
         if (ds<=dd)&&(y0<=y1_coor+hh/2)&&(y0>=y1_coor-hh/2)
          nodesc = [nodesc;co2id([x,y,z])]; 
         end
      end
   end
end

case 'z'
for z = 1:kz+1
   for y = 1:ky+1
      for x = 1:kx+1
         tmpcounter = 1 + tmpcounter;
	 x0=nodes(tmpcounter,1);
	 y0=nodes(tmpcounter,2);
	 z0=nodes(tmpcounter,3);
         ds=sqrt((x0-x1_coor)^2+(y0-y1_coor)^2);
         if (ds<=dd)&&(z0<=z1_coor+hh/2)&&(z0>=z1_coor-hh/2)
          nodesc = [nodesc;co2id([x,y,z])]; 
         end
      end
   end
end

end
