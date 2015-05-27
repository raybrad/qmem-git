function [currentT,currddT] = tdcurrentf(cod,dir,nsteps,qmcuurr)

global sigma epsilon_in epsilon_sd;
global epsilon_mt;
global mun mup;
global scl;
global nodes links contacts;
global Nnode Nlink;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
	
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes dirNodes;
global isBndNodes isDirNodes dcVolDirNodes; %#ok<*NUSED>
global doping; % doping profile
global semiNodes eqnSemiNodes; %#ok<NUSED>
global dirSemiNodes isDirSemiNodes;
global kx ky kz;
global epsilon_mt;
global savefile;
global isqmvolm;
global epsilon_qm;

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
currentT = zeros(nsteps,1);
currddT  =zeros(length(Snodes(:,1)),(nsteps+3));

for TT=1:nsteps

filename = [savefile,num2str(TT-1),'.mat'];
load(filename);

 for i = 1:length(Snodes(:,1))
	n1=Snodes(i,1);
	n2=Snodes(i,2);

        currddT(i,1) = (nodes(n1,1)+nodes(n2,1))/2;
        currddT(i,2) = (nodes(n1,2)+nodes(n2,2))/2;
        currddT(i,3) = (nodes(n1,3)+nodes(n2,3))/2;

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
        for qmi = 1:length(slkvol)
          if isqmvolm(slkvol(qmi))
           slkvolm(qmi)=4;
          end
        end
	slinkl  = linkL(slink);
	alpha_n = mun/slinkl; alpha_p = mup/slinkl;

        dtE1 = -(dtVp(n2)-dtVp(n1))/slinkl - dtHp(slink);
   
        if isDirSemiNodes(n1), dV1 = deltaV(doping(n1)); else dV1 = 0; end
        if isDirSemiNodes(n2), dV2 = deltaV(doping(n2)); else dV2 = 0; end
        
        if all(slkvolm == 4)

          J = qmcuurr(TT) + epsilon_qm *dtE1;
         
          currddT(i,(TT+3)) = currddT(i,(TT+3))+ J;
          currentT(TT) = currentT(TT) + linkS(slink) * J;

        else 
	 for t = 1:length(slkvol)
		switch slkvolm(t)
			case 1
                                beta = V(n2)-V(n1)+H(slink)*slinkl;
                                J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                                       +Jc(alpha_p,beta,'p',p(n1),p(n2));
                                J = J_diff + epsilon_sd*dtE1;
			case 2
                                J = sigma*((V(n1)-V(n2)+dV1-dV2)/slinkl - H(slink))+ epsilon_mt*dtE1;
			case 3
                                J = epsilon_in*dtE1;
                        case 4
                                beta = V(n2)-V(n1)+H(slink)*slinkl;
                                J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                                       +Jc(alpha_p,beta,'p',p(n1),p(n2));
                                J = J_diff + epsilon_qm*dtE1;
			otherwise
			   error('undefined material');
		end	
    
          currddT(i,(TT+3)) = currddT(i,(TT+3))+ J; 
          currentT(TT) = currentT(TT) + slkvols(t) * J;
         end
        end
 end
end
