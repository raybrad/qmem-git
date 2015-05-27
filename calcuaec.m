function abe = calcuae(V,H)

global sigma epsilon_in epsilon_sd;
global epsilon_mt;
global mun mup;
global scl;
global nodes links contacts;
global Nnode Nlink Nsurf;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes dirNodes;
global isBndNodes isDirNodes dcVolDirNodes acVolDirNodes;
global doping; 
global semiNodes;
global dirSemiNodes isDirSemiNodes;
global kx ky kz;
global bndLinks;
global metalNodes;
global savefile;
global isSQMlinks;

%%%%%%% initial guess  %%%%%%%%%%%%%
EL  = zeros(Nlink,1);

E = zeros(Nnode,4);

eqnLinks = (1:Nlink)';
eqnNodes = (1:Nnode)';
eqnSurfs = (1:Nsurf)';

%%%%%%%% start Newton's iteration %%%%%%%%
newtonTol = 1e-15;
updateTol = 1e-10;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 20;
maxLinIt = 200;
Eps = [epsilon_sd,epsilon_mt,epsilon_in];
Sgm = [0,sigma,0];

Nl = length(eqnLinks);
NeqnNodes = length(eqnNodes);
Ns = length(eqnSurfs);

normRes = 1;
normUpdate = 1;
normRes_pre = 1e10;
itNr = 0;

    %%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
     for l1 = eqnLinks'
      if isSQMlinks(l1) == 0
        n1 = links(l1,1);
        n2 = links(l1,2);
        n3 = links(l1,3);
        if nodes(n1,n3) < nodes(n2,n3)
          n1f = n1;
          n2f = n2;
        else
          n1f = n2;
          n2f = n1;
        end
        ajvol_l1 = linkVolumes{l1}(1,:);
        ajvolS_l1 = linkVolumes{l1}(2,:);
        ajvolM_l1 = volumeM(ajvol_l1);

        dV_vec = [0,1,1/2];

        if isDirSemiNodes(n1f), dV1 = deltaV(doping(n1f)); else dV1 = 0; end
        if isDirSemiNodes(n2f), dV2 = deltaV(doping(n2f)); else dV2 = 0; end

        ETT=sum(Eps(ajvolM_l1).*(-(V(n2f)-V(n1f)+dV_vec(ajvolM_l1).*dV2-dV_vec(ajvolM_l1).*dV1)/linkL(l1)-H(l1)).*ajvolS_l1);
        AAT=sum(Eps(ajvolM_l1).*ajvolS_l1);
        if AAT ~= 0
         EL(l1) = ETT/AAT;
        end
      end
     end

     for n1 = eqnNodes'
         ett = zeros(3,2);
         ajlk_n1 = nodeLinks{n1}(1,:);

         for i = 1:length(ajlk_n1)
             lk = ajlk_n1(i);           
             ett(links(lk,3),1) = ett(links(lk,3),1)+ linkL(lk);
             ett(links(lk,3),2) = ett(links(lk,3),2)+ EL(lk)*linkL(lk);             
         end
        if ett(1,1) ~=0
           E(n1,1) = ett(1,2)/ett(1,1);
        end
        if ett(2,1) ~=0
           E(n1,2) = ett(2,2)/ett(2,1);
        end
        if ett(3,1) ~=0
           E(n1,3) = ett(3,2)/ett(3,1);
        end
           
        E(n1,4)=sqrt(E(n1,1)^2 + E(n1,2)^2 + E(n1,3)^2);

     end

     abe = E(:,4);
    

