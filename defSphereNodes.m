function nodesc = defSphereNodes(p1,lambda,radius,nodes)

%%define the nodes of a sphere
%%p1: the center coordinate of sphere
%%radius: the radius 

global kx ky kz;
x1 = p1(1)*1e-9/lambda; y1 = p1(2)*1e-9/lambda; z1 = p1(3)*1e-9/lambda;
radius=radius*1e-9/lambda;


nodesc = [];

tmpcounter = 0;
for z = 1:kz+1
   for y = 1:ky+1
      for x = 1:kx+1
         tmpcounter = 1 + tmpcounter;
	 x0=nodes(tmpcounter,1);
	 y0=nodes(tmpcounter,2);
	 z0=nodes(tmpcounter,3);
         ds=sqrt((x0-x1)^2+(y0-y1)^2+(z0-z1)^2);
         if ds<=radius
          nodesc = [nodesc;co2id([x,y,z])]; 
         end
      end
   end
end


