function [currentdqm,currentdqmn] = currentSqm(Vs,ns,ps)

global sigma epsilon_in epsilon_sd;
global mun mup;
global scl;
global nodes links contacts;
global Nnode Nlink;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes dirNodes;
global isBndNodes isDirNodes dcVolDirNodes; 
global doping; % doping profile
global semiNodes; 
global dirSemiNodes isDirSemiNodes;
global kx ky kz;
global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global sQMlinks;
global isqmnodes;

global ehcrf;


currentdqm  = zeros(Nlink,6);
currentdqmn = zeros(Nnode,6);

for ilk =1:Nlink
    n1  = links(ilk,1); n2 = links(ilk,2); n3 = links(ilk,3);
    currentdqm(ilk,1)=(nodes(n1,1)+nodes(n2,1))/2;
    currentdqm(ilk,2)=(nodes(n1,2)+nodes(n2,2))/2;
    currentdqm(ilk,3)=(nodes(n1,3)+nodes(n2,3))/2;
    if  isSQMlinks(ilk) == 0;

        v1  = Vs(n1);
        v2  = Vs(n2);
        ns1 = ns(n1);
        ns2 = ns(n2);
        ps1 = ps(n1);
        ps2 = ps(n2);
    
        if isDirSemiNodes(n1)  dV1 = deltaV(doping(n1)); else dV1 = 0; end
        if isDirSemiNodes(n2)  dV2 = deltaV(doping(n2)); else dV2 = 0; end

        slkvol  = linkVolumes{ilk}(1,:);
        slkvols = linkVolumes{ilk}(2,:);
        slkvolm = volumeM(slkvol);
        slinkl  = linkL(ilk);
        alpha_n = mun/slinkl; alpha_p = mup/slinkl;

        for t = 1:length(slkvol)
                switch slkvolm(t)
                        case 1
                             beta = v2-v1;
                             J = Jc(alpha_n,beta,'n',ns1,ns2)...
                                 +Jc(alpha_p,beta,'p',ps1,ps2);
                        case 2
                             J = sigma*(v1-v2 + dV1-dV2)/slinkl;
                        case 3
                             J = 0;
                        otherwise
                           error('undefined material');
                end
          currentdqm(ilk,n3+3) = currentdqm(ilk,n3+3)+ J;
        end
    end
end

     for n1 = 1:Nnode

    currentdqmn(n1,1)=nodes(n1,1);
    currentdqmn(n1,2)=nodes(n1,2);
    currentdqmn(n1,3)=nodes(n1,3);
          
         ett = zeros(3,2);
         ajlk_n1 = nodeLinks{n1}(1,:);
          
         for i = 1:length(ajlk_n1)
             lk = ajlk_n1(i);
             ett(links(lk,3),1) = ett(links(lk,3),1)+ linkL(lk);
             ett(links(lk,3),2) = ett(links(lk,3),2)+ currentdqm(lk,links(lk,3)+3)*linkL(lk);
         end
        if ett(1,1) ~=0
           currentdqmn(n1,4) = ett(1,2)/ett(1,1);
        end
        if ett(2,1) ~=0
           currentdqmn(n1,5) = ett(2,2)/ett(2,1);
        end
        if ett(3,1) ~=0
           currentdqmn(n1,6) = ett(3,2)/ett(3,1);
        end
end

