function plotStructure(metalNodes,Nnode,kx,ky,kz,x_coor,y_coor,z_coor)

ismetalNodes = zeros(Nnode,1);
ismetalNodes(metalNodes)=1;
tmpcounter=0;
isTranmetalNodes(1:kx+1,1:ky+1,1:kz+1)=0;
for k = 1:kz+1
    for j = 1:ky+1
        for i = 1:kx+1
	tmpcounter=tmpcounter+1
	isTranmetalNodes(i,j,k)=ismetalNodes(tmpcounter);
	end
    end
end
[X,Y,Z]=ndgrid(x_coor,y_coor,z_coor);

figure;
scatter3(X,Y,Z,isTranmetalNodes);
exit;
