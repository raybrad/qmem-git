function nodes = defBrickNodes(p1,p2)

%%% provide all the nodes within a cube defined by its lower left and upper right
%%% points

x1 = p1(1); y1 = p1(2); z1 = p1(3);
x2 = p2(1); y2 = p2(2); z2 = p2(3);

nodes = [];
for z = z1:z2
   for y = y1:y2
      for x = x1:x2
         nodes = [nodes;co2id([x,y,z])]; 
      end
   end
end



% if xn == 0
%     ynodes = (id1:(kx+1):id1+yn*(kx+1))';
%     for i = 0:zn
%         nodes = [nodes;ynodes+i*(kx+1)*(ky+1)];
%     end
% elseif yn == 0
%     xnodes = (id1:(id1+xn))';
%     for i = 0:zn
%        nodes = [nodes;xnodes+i*(kx+1)*(ky+1)];
%     end
% elseif zn == 0
%     xnodes = (id1:(id1+xn))';
%     for i = 0:yn
%        nodes = [nodes;xnodes+i*(kx+1)];
%     end
% else
%     error('Non-Cartesian brisk');
% end