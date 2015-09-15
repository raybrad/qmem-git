function outputField(step,outputPosCom,outputPlane,V,A,H)

global scl;
global dt;
global x_coor y_coor z_coor;
global Origin Radius;

no_of_nodes_x = length(x_coor);   %number of nodes along x-axis
no_of_nodes_y = length(y_coor);   %number of nodes along y-axis
no_of_nodes_z = length(z_coor);   %number of nodes along z-axis

lunitx=x_coor*1e-9/scl.lambda;
lunity=y_coor*1e-9/scl.lambda;
lunitz=z_coor*1e-9/scl.lambda;

DiVecPotx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
DiVecPoty(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
DiVecPotz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;

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
             	 DiVecPotx(i,j,k)=H(tmpcounter);
             	 VecPotx(i,j,k)=A(tmpcounter);
		        else
	           	 DiVecPotx(i,j,k)=DiVecPotx(i-1,j,k);
	           	 VecPotx(i,j,k)=VecPotx(i-1,j,k);
	        	end
		        if (i>1 && i<no_of_nodes_x)
		         DiVecPotx(i,j,k)=0.5*(DiVecPotx(i,j,k)+DiVecPotx(i-1,j,k));
		         VecPotx(i,j,k)=0.5*(VecPotx(i,j,k)+VecPotx(i-1,j,k));
		        end

                if (j < no_of_nodes_y) 
                 tmpcounter = 1 + tmpcounter;
		         DiVecPoty(i,j,k)=H(tmpcounter);
		         VecPoty(i,j,k)=A(tmpcounter);
	        	else
	          	 DiVecPoty(i,j,k)=DiVecPoty(i,j-1,k);
	          	 VecPoty(i,j,k)=VecPoty(i,j-1,k);
		        end
     		    if (j>1 && j<no_of_nodes_y)
		         DiVecPoty(i,j,k)=0.5*(DiVecPoty(i,j,k)+DiVecPoty(i,j-1,k));
		         VecPoty(i,j,k)=0.5*(VecPoty(i,j,k)+VecPoty(i,j-1,k));
		        end

                if (k < no_of_nodes_z) 
                 tmpcounter = 1 + tmpcounter;
		         DiVecPotz(i,j,k)=H(tmpcounter);
		         VecPotz(i,j,k)=A(tmpcounter);
		        else
		         DiVecPotz(i,j,k)=DiVecPotz(i,j,k-1);
		         VecPotz(i,j,k)=VecPotz(i,j,k-1);
		        end
		        if (k>1 && k<no_of_nodes_z)
		         DiVecPotz(i,j,k)=0.5*(DiVecPotz(i,j,k)+DiVecPotz(i,j,k-1));
		         VecPotz(i,j,k)=0.5*(VecPotz(i,j,k)+VecPotz(i,j,k-1));
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
[dVy,dVx,dVz]=gradient(Voltage,lunity,lunitx,lunitz);

Ex(:,:,:)=-dVx(:,:,:)*scl.s_E-DiVecPotx(:,:,:)*scl.s_E;
Ey(:,:,:)=-dVy(:,:,:)*scl.s_E-DiVecPoty(:,:,:)*scl.s_E;
Ez(:,:,:)=-dVz(:,:,:)*scl.s_E-DiVecPotz(:,:,:)*scl.s_E;

[py,px,pz]=gradient(VecPotx,lunity,lunitx,lunitz);
[qy,qx,qz]=gradient(VecPoty,lunity,lunitx,lunitz);
[ry,rx,rz]=gradient(VecPotz,lunity,lunitx,lunitz);

Bx = ry-qz;
By = pz-rx;
Bz = qx-py;

Bx=Bx*scl.s_B;
By=By*scl.s_B;
Bz=Bz*scl.s_B;

if (step==1) 
    delete('*Ex.dat','*E2.dat');
    delete('*By.dat','*B2.dat');
end
%output point
for i=1:length(outputPosCom)
    
    filename=outputPosCom{i,3};
    %fprintf('filename:%s \n',filename);
    fpID=fopen(filename,'a');
    idx=outputPosCom{i,1}(1);
    idy=outputPosCom{i,1}(2);
    idz=outputPosCom{i,1}(3);
    switch outputPosCom{i,2}
    case 'Ex'
        out_point=Ex(idx,idy,idz);
    case 'Ey'
        out_point=Ey(idx,idy,idz);
    case 'Ez'
        out_point=Ez(idx,idy,idz);
    case 'Bx'
        out_point=Bx(idx,idy,idz);
    case 'By'
        out_point=By(idx,idy,idz);
    case 'Bz'
        out_point=Bz(idx,idy,idz);
    otherwise
        display('error,component not fould');
    end
        %fprintf('E %12.5e %12.5e \n',dt*scl.tao*step/1e-15,out_point);
        fprintf(fpID,'%12.5e %12.5e \n',dt*scl.tao*step/1e-15,out_point);
        fclose(fpID);
end

%output plane
if (mod(step,100)==0) 
out_plane(1:no_of_nodes_x,1:no_of_nodes_y) =0.0;
filename=strcat(num2str(step),outputPlane{1,3});
fpID=fopen(filename,'w');
idp=outputPlane{1,2};
switch outputPlane{1,1}
case 'z'
    out_plane=sqrt(Ex(:,:,idp).^2+Ey(:,:,idp).^2);
case 'y'
    out_plane=sqrt(Ex(:,idp,:).^2+Ez(:,idp,:).^2);
case 'x'
    out_plane=sqrt(Ey(idp,:,:).^2+Ez(idp,:,:).^2);
end

for i=1:no_of_nodes_x
for j=1:no_of_nodes_y
fprintf(fpID,'%12.5e', out_plane(i,j));
end
fprintf(fpID,'\n');
end
fclose(fpID);
end

