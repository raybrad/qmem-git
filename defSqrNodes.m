function nodesc = defSqrNodes(p1,p2,sideLength,height,ll,ww,nodes,lambda)

%%define the nodes of a bowtie antenna x directed
%%height: the thickness
global kx ky kz;
x1 = p1(1)*1e-9/lambda; y1 = p1(2)*1e-9/lambda; z1 = p1(3)*1e-9/lambda;
x2 = p2(1)*1e-9/lambda; y2 = p2(2)*1e-9/lambda; z2 = p2(3)*1e-9/lambda;
sideLength=sideLength*1e-9/lambda;
height=height*1e-9/lambda;
ll=ll*1e-9/lambda;
ww=ww*1e-9/lambda;

nodesc = [];

tmpcounter = 0;
for z = 1:kz+1
   for y = 1:ky+1
      for x = 1:kx+1
         tmpcounter = 1 + tmpcounter;
	 x0=nodes(tmpcounter,1);
	 y0=nodes(tmpcounter,2);
	 z0=nodes(tmpcounter,3);
%big block part
	  if( y0>=y1-sideLength/2.0 && y0<=y1+sideLength/2.0 && z0>=z1-height/2.0 && z0<z1+height/2.0)
		if ((x0>=x1 && x0<=x2) || (x0>=x2 && x0<=x1))
             nodesc = [nodesc;co2id([x,y,z])]; 
	  end
	end
%small block part
	  if( y0>= y2-ww/2.0 && y0<= y2+ ww/2.0 && z0>= z2-ww/2.0 && z0<= z2+ ww/2.0)
	     if (x2>x1  && x0>=x2 && x0 <=x2+ll)
             nodesc = [nodesc;co2id([x,y,z])]; 
	     elseif (x2<x1  && x0<=x2 && x0 >=x2-ll)
             nodesc = [nodesc;co2id([x,y,z])]; 
	     end
	  end
      end
   end
end



