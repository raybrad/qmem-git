function poyvector(axisxyz,m)
fname=strcat('dump/variables_',num2str(m),'.mat');
load(fname);
 
no_of_nodes_x=16;
no_of_nodes_y=16;
no_of_nodes_z=16;

%read dA/dt
DiVecPotx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
DiVecPoty(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
DiVecPotz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;

tmpcounter = 0;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
                if (i < no_of_nodes_x) 
	        tmpcounter = 1 + tmpcounter;
		DiVecPotx(i,j,k)=H(tmpcounter);
		else
		DiVecPotx(i,j,k)=DiVecPotx(i-1,j,k);
		end
		   if (i>1 && i<no_of_nodes_x)
		   DiVecPotx(i,j,k)=0.5*(DiVecPotx(i,j,k)+DiVecPotx(i-1,j,k));
		   end

                if (j < no_of_nodes_y) 
                tmpcounter = 1 + tmpcounter;
		DiVecPoty(i,j,k)=H(tmpcounter);
		else
		DiVecPoty(i,j,k)=DiVecPoty(i,j-1,k);
		end
		   if (j>1 && j<no_of_nodes_y)
		   DiVecPoty(i,j,k)=0.5*(DiVecPoty(i,j,k)+DiVecPoty(i,j-1,k));
		   end

                if (k < no_of_nodes_z) 
                tmpcounter = 1 + tmpcounter;
		DiVecPotz(i,j,k)=H(tmpcounter);
		else
		DiVecPotz(i,j,k)=DiVecPotz(i,j,k-1);
		end
		   if (k>1 && k<no_of_nodes_z)
		   DiVecPotz(i,j,k)=0.5*(DiVecPotz(i,j,k)+DiVecPotz(i,j,k-1));
		   end
        end
    end
end
% get V and and div V
tmpcounter = 0;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
            tmpcounter = (no_of_nodes_x*(j-1))+(no_of_nodes_y*no_of_nodes_x*(k-1))+i;
            Voltage(i,j,k)=V(tmpcounter);
            
        end
    end
end
% column first , then rows
[dVy,dVx,dVz]=gradient(Voltage,1.0,1.0,1.0);

Ex(:,:,:)=-dVx(:,:,:)-DiVecPotx(:,:,:);
Ey(:,:,:)=-dVy(:,:,:)-DiVecPoty(:,:,:);
Ez(:,:,:)=-dVz(:,:,:)-DiVecPotz(:,:,:);

% get B field
VecPotx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
VecPoty(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
VecPotz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;

Bx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
By(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
Bz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;

tmpcounter = 0;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
                if (i < no_of_nodes_x) 
	        tmpcounter = 1 + tmpcounter;
		VecPotx(i,j,k)=A(tmpcounter);
		else
		VecPotx(i,j,k)=VecPotx(i-1,j,k);
		end
		   if (i>1 && i<no_of_nodes_x)
		   VecPotx(i,j,k)=0.5*(VecPotx(i,j,k)+VecPotx(i-1,j,k));
		   end

                if (j < no_of_nodes_y) 
                tmpcounter = 1 + tmpcounter;
		VecPoty(i,j,k)=A(tmpcounter);
		else
		VecPoty(i,j,k)=VecPoty(i,j-1,k);
		end
		   if (j>1 && j<no_of_nodes_y)
		   VecPoty(i,j,k)=0.5*(VecPoty(i,j,k)+VecPoty(i,j-1,k));
		   end

                if (k < no_of_nodes_z) 
                tmpcounter = 1 + tmpcounter;
		VecPotz(i,j,k)=A(tmpcounter);
		else
		VecPotz(i,j,k)=VecPotz(i,j,k-1);
		end
		   if (k>1 && k<no_of_nodes_z)
		   VecPotz(i,j,k)=0.5*(VecPotz(i,j,k)+VecPotz(i,j,k-1));
		   end
        end
    end
end
%
%[Bx,By,Bz]=curl(VecPotx,VecPoty,VecPotz);
[py,px,pz]=gradient(VecPotx);
[qy,qx,qz]=gradient(VecPoty);
[ry,rx,rz]=gradient(VecPotz);

Bx = ry-qz;
By = pz-rx;
Bz = qx-py;


%get poynting vector
poyvecx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
poyvecy(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
poyvecz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
	    poyvecx(i,j,k)=  Ey(i,j,k)*Bz(i,j,k)-Ez(i,j,k)*By(i,j,k);
	    poyvecy(i,j,k)= -Ex(i,j,k)*Bz(i,j,k)+Ez(i,j,k)*Bx(i,j,k);
	    poyvecz(i,j,k)=  Ex(i,j,k)*By(i,j,k)-Ey(i,j,k)*Bx(i,j,k);
	end
    end
end




%[X,Y]=ndgrid(1.0:1.0:16.0,1.0:1.0:16.0);

%switch axisxyz
%case 1
%EE1(:,:)=poyvecy(8,:,:);
%EE2(:,:)=poyvecz(8,:,:);
%case 2
%EE1(:,:)=poyvecx(:,8,:);
%EE2(:,:)=poyvecz(:,8,:);
%case 3 
%EE1(:,:)=poyvecx(:,:,8);
%EE2(:,:)=poyvecy(:,:,8);
%otherwise
%display('error');
%end


%%quiver 3D vector arrows
%figure;
%quiver(X,Y,EE1,EE2);


figure;
[X,Y,Z]=ndgrid(1.0:1.0:16.0,1.0:1.0:16.0,1.0:1.0:16.0);
Exx(:,:,:)=poyvecx(1:1:16,1:1:16,1:1:16);
Eyy(:,:,:)=poyvecy(1:1:16,1:1:16,1:1:16);
Ezz(:,:,:)=poyvecz(1:1:16,1:1:16,1:1:16);

quiver3(X,Y,Z,Exx,Eyy,Ezz,'color',[0,0,1]);

clear;
