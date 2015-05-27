function [current,currentd,totsarea] = avcurrent(cod,Vs,ns,ps,V,n,p,A,dir)

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

for i = 1:length(Snodes(:,1))
	n1=Snodes(i,1);
	n2=Snodes(i,2);
	adlk = nodeLinks{n1}(1,:);
	adnd = nodeLinks{n1}(2,:);
	for j= 1:length(adlk)
		if adnd(j) == n2
			slink = adlk(j);
		end
	end
    
	Vs1  = Vs(n1);
	Vs2  = Vs(n2);
	ns1  = ns(n1);
	ns2  = ns(n2);
	ps1  = ps(n1);
	ps2  = ps(n2);

        slkvol  = linkVolumes{slink}(1,:);
	slkvols = linkVolumes{slink}(2,:);
	slkvolm = volumeM(slkvol);
	slinkl  = linkL(slink);
	alpha_n = mun/slinkl; alpha_p = mup/slinkl;
   
        if isDirSemiNodes(n1), v1 = 0; else v1 = 1; end
        if isDirSemiNodes(n2), v2 = 0; else v2 = 1; end
        beta = Vs2-Vs1;
        E1 = -(V(n2)-V(n1))/slinkl-1i*omega*A(slink);

	for t = 1:length(slkvol)
		switch slkvolm(t)
			case 1
                                sigma_n = sigma_c(mun,beta,'n',ns1,ns2);
                                sigma_p = sigma_c(mup,beta,'p',ps1,ps2);

                                J_diff = Jc1_diff(alpha_n,beta,'n',n(n1),n(n2),v1,v2)...
                                    +Jc1_diff(alpha_p,beta,'p',p(n1),p(n2),v2,v2);
                                J_drift = (sigma_n+sigma_p)*E1;
                                J = J_diff+J_drift;
			case 2
                                J = sigma*E1;
			case 3
                                J = 0;
			otherwise
			   error('undefined material');
		end	
     
        current = current + slkvols(t) * J;
        totsarea= totsarea + slkvols(t);
    end
end

if (totsarea ~= 0)
	currentd = current/totsarea;
end

