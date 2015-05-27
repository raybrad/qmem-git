function nodesc = defTriNodes(p1,p2,p3,bottom,top,height,ll,ww,nodes,lambda)

%%define the nodes of a bowtie antenna x directed
%%p1: the center coordinate of bottom line
%%p2: the center coordinate of topi   line
%%height: the thickness
global kx ky kz;
x1 = p1(1)*1e-9/lambda; y1 = p1(2)*1e-9/lambda; z1 = p1(3)*1e-9/lambda;
x2 = p2(1)*1e-9/lambda; y2 = p2(2)*1e-9/lambda; z2 = p2(3)*1e-9/lambda;
x3 = p3(1)*1e-9/lambda; y3 = p3(2)*1e-9/lambda; z3 = p3(3)*1e-9/lambda;
bottom=bottom*1e-9/lambda;
top=top*1e-9/lambda;
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
%Tri part
     	 if ( z0>=z1-height/2.0 && z0<=z1+height/2.0)
	  if ( x0>= x1 &&  x0<= x2 )
	    temp= (x2-x0)/(x2-x1)*(bottom-top)+top;
	    if (y0>=y1-temp/2.0 && y0<=y1+temp/2.0)
             nodesc = [nodesc;co2id([x,y,z])]; 
	     end
	  elseif(x0>=x2 && x0<= x1)
	     temp= (x0-x2)/(x1-x2)*(bottom-top)+top;
	     if (y0>=y1-temp/2.0 && y0<=y1+temp/2.0)
             nodesc = [nodesc;co2id([x,y,z])]; 
	     end
	  end
        end
%big block part
	  if( y0>=y2-top/2.0 && y0<=y2+top/2.0 && z0>=z2-height/2.0 && z0<z2+height/2.0)
		if ((x0>=x2 && x0<=x3) || (x0>=x3 && x0<=x2))
             nodesc = [nodesc;co2id([x,y,z])]; 
	  end
	end
%small block part
	  if( y0>= y3-ww/2.0 && y0<= y3+ ww/2.0 && z0>= z3-ww/2.0 && z0<= z3+ ww/2.0)
	     if (x3>x1  && x0>=x3 && x0 <=x3+ll)
             nodesc = [nodesc;co2id([x,y,z])]; 
	     elseif (x3<x1  && x0<=x3 && x0 >=x3-ll)
             nodesc = [nodesc;co2id([x,y,z])]; 
	     end
	  end
      end
   end
end



