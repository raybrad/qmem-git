function [current,currentd,totsarea,currdd] = currentStd(cod,V,n,p,H,dtV,dtH,dir)

global sigma epsilon_in epsilon_sd;
global mun mup;
global scl;
global nodes links contacts;
global Nnode;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
	
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes dirNodes;
global isBndNodes isDirNodes dcVolDirNodes; %#ok<*NUSED>
global doping; % doping profile
global semiNodes eqnSemiNodes; %#ok<NUSED>
global dirSemiNodes isDirSemiNodes;
global kx ky kz;
global omega epsilon_mt;
global savefile;

Snodes1 = [];
Snodes2 = [];

if (dir == 1)
	nx  = cod(1);
	for z = 1 : kz+1
		for y = 1 : ky+1
			Snodes1 = [Snodes1;co2id([nx,y,z])];
			Snodes2 = [Snodes2;co2id([nx+1,y,z])];
		end
	end
elseif (dir == 2)
    ny  = cod(2);
    for z = 1 : kz+1
		for x = 1 : kx+1
		    Snodes1 = [Snodes1;co2id([x,ny,z])];
			Snodes2 = [Snodes2;co2id([x,ny+1,z])];
		end
	end
elseif (dir == 3)
    nz  = cod(3);
	for y = 1 : ky+1
		for x = 1 : kx+1
			Snodes1 = [Snodes1;co2id([x,y,nz])];
			Snodes2 = [Snodes2;co2id([x,y,nz+1])];
		end
	end
end

Snodes  = [Snodes1,Snodes2];
current = 0;
currentd= 0;
totsarea= 0;
currdd  =zeros(length(Snodes(:,1)),4);

for i = 1:length(Snodes(:,1))
	n1=Snodes(i,1);
	n2=Snodes(i,2);

        currdd(i,1) = (nodes(n1,1)+nodes(n2,1))/2;
        currdd(i,2) = (nodes(n1,2)+nodes(n2,2))/2;
        currdd(i,3) = (nodes(n1,3)+nodes(n2,3))/2;

	adlk = nodeLinks{n1}(1,:);
	adnd = nodeLinks{n1}(2,:);
	for j= 1:length(adlk)
		if adnd(j) == n2
			slink = adlk(j);
		end
	end
    

        slkvol  = linkVolumes{slink}(1,:);
	slkvols = linkVolumes{slink}(2,:);
	slkvolm = volumeM(slkvol);
	slinkl  = linkL(slink);
	alpha_n = mun/slinkl; alpha_p = mup/slinkl;
   
        if isDirSemiNodes(n1), dV1 = deltaV(doping(n1)); else dV1 = 0; end
        if isDirSemiNodes(n2), dV2 = deltaV(doping(n2)); else dV2 = 0; end

        dtE1 = -(dtV(n2)-dtV(n1))/slinkl - dtH(slink);

	for t = 1:length(slkvol)
		switch slkvolm(t)
			case 1
                                beta = V(n2)-V(n1)+H(slink)*slinkl;
                                J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                                       +Jc(alpha_p,beta,'p',p(n1),p(n2));
                                J = J_diff + epsilon_sd*dtE1;
			case 2
                                J = sigma*((V(n1)-V(n2)+dV1-dV2)/slinkl - H(slink)) + epsilon_mt*dtE1;
			case 3
                                J = epsilon_in*dtE1;
			otherwise
			   error('undefined material');
		end	
    
        currdd(i,4) = currdd(i,4)+ J; 
        current = current + slkvols(t) * J;
        totsarea= totsarea + slkvols(t);
    end
end

if (totsarea ~= 0)
	currentd = current/totsarea;
end

