function [E,B,AN] = tdcalcueb(V,n,p,A,H)

display(['  Start TD calculation for E and B:']);

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

%%%%%%% initial guess  %%%%%%%%%%%%%
EL  = zeros(Nlink,1);
BS  = zeros(Nsurf,2);

E = zeros(Nnode,4);
B = zeros(Nnode,4);
AN= zeros(Nnode,4);

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
tic;

    %%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
     for l1 = eqnLinks'
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

     for s1 = eqnSurfs' 
         l1 = surfLinks(s1,1); 
         l2 = surfLinks(s1,2); 
         l3 = surfLinks(s1,3); 
         l4 = surfLinks(s1,4);

         sare = linkL(l1)*linkL(l2);

         n11 = links(l1,1);
         n12 = links(l1,2);

         n21 = links(l2,1);
         n22 = links(l2,2);
          
         n31 = links(l1,3);
         n32 = links(l2,3);

        if n31==1 && n32 ==3
          if nodes(n11,n31) < nodes(n12,n31) && nodes(n21,n32) < nodes(n22,n32)
           BS(s1,1)=-(A(l1)*linkL(l1)+A(l2)*linkL(l2)-A(l3)*linkL(l3)-A(l4)*linkL(l4))/sare;
          else 
           error('circule1');
          end
        elseif n31==1 && n32 ==2
          if nodes(n11,n31) < nodes(n12,n31) && nodes(n21,n32) < nodes(n22,n32)
           BS(s1,1)=(A(l1)*linkL(l1)+A(l2)*linkL(l2)-A(l3)*linkL(l3)-A(l4)*linkL(l4))/sare;
          else
           error('circule2');
          end
        elseif n31==2 && n32 ==3
          if nodes(n11,n31) < nodes(n12,n31) && nodes(n21,n32) < nodes(n22,n32)
           BS(s1,1)=(A(l1)*linkL(l1)+A(l2)*linkL(l2)-A(l3)*linkL(l3)-A(l4)*linkL(l4))/sare;
          else
           error('circule3');
          end
        else 
           error('circule4');             
        end

         if (n31 ~= 1) && (n32 ~= 1)
           BS(s1,2)=1;
         elseif (n31 ~= 2) && (n32 ~= 2)
           BS(s1,2)=2;
         elseif (n31 ~= 3) && (n32 ~= 3)
           BS(s1,2)=3;
         end
        
     end

     for n1 = eqnNodes'
         
         ett = zeros(3,3);
         btt = zeros(3,2); 
         ajlk_n1 = nodeLinks{n1}(1,:);

         for i = 1:length(ajlk_n1)
             lk = ajlk_n1(i);           
             ett(links(lk,3),1) = ett(links(lk,3),1)+ linkL(lk);
             ett(links(lk,3),2) = ett(links(lk,3),2)+ EL(lk)*linkL(lk);             
             ett(links(lk,3),3) = ett(links(lk,3),3)+ A(lk)*linkL(lk);
             ajsf_l1 = linkSurfs{lk};  
            
             for ii = 1:size(ajsf_l1,2)
              ajlk = ajsf_l1(2:4,ii);
              if all(links(ajlk,3) ~= 1)
               if BS(ajsf_l1(1,ii),2)~=1
                 error('surface v1');
               end
               btt(1,1)= btt(1,1)+ linkL(ajlk(1))*linkL(ajlk(2));
               btt(1,2)= btt(1,2)+BS(ajsf_l1(1,ii),1)*linkL(ajlk(1))*linkL(ajlk(2));
              elseif all(links(ajlk,3) ~= 2)
               if BS(ajsf_l1(1,ii),2)~=2
                 error('surface v2');
               end
               btt(2,1)= btt(2,1)+linkL(ajlk(1))*linkL(ajlk(2));
               btt(2,2)= btt(2,2)+BS(ajsf_l1(1,ii),1)*linkL(ajlk(1))*linkL(ajlk(2));
              elseif all(links(ajlk,3) ~= 3)
               if BS(ajsf_l1(1,ii),2)~=3
                 error('surface v3');
               end
               btt(3,1)= btt(3,1)+linkL(ajlk(1))*linkL(ajlk(2));
               btt(3,2)= btt(3,2)+BS(ajsf_l1(1,ii),1)*linkL(ajlk(1))*linkL(ajlk(2));
              end
             end
         end
        if ett(1,1) ~=0

           E(n1,1) = ett(1,2)/ett(1,1);
           AN(n1,1) = ett(1,3)/ett(1,1);
        end
        if ett(2,1) ~=0

           E(n1,2) = ett(2,2)/ett(2,1);
           AN(n1,2) = ett(2,3)/ett(2,1);
        end
        if ett(3,1) ~=0

           E(n1,3) = ett(3,2)/ett(3,1);
           AN(n1,3) = ett(3,3)/ett(3,1);
        end
           
        E(n1,4)=sqrt(E(n1,1)^2 + E(n1,2)^2 + E(n1,3)^2);
        AN(n1,4)=sqrt(AN(n1,1)^2 + AN(n1,2)^2 + AN(n1,3)^2); 
        if E(n1,4) ~= 0
           E(n1,1)=E(n1,1)/E(n1,4);
           E(n1,2)=E(n1,2)/E(n1,4);
           E(n1,3)=E(n1,3)/E(n1,4);
        end

        if AN(n1,4) ~= 0
           AN(n1,1)=AN(n1,1)/AN(n1,4);
           AN(n1,2)=AN(n1,2)/AN(n1,4);
           AN(n1,3)=AN(n1,3)/AN(n1,4);
        end

        if btt(1,1) ~=0

           B(n1,1) = btt(1,2)/btt(1,1);
        end
        if btt(2,1) ~=0

           B(n1,2) = btt(2,2)/btt(2,1);
        end
        if btt(3,1) ~=0

           B(n1,3) = btt(3,2)/btt(3,1);
        end

        B(n1,4)=sqrt(B(n1,1)^2 + B(n1,2)^2 + B(n1,3)^2);

        if B(n1,4) ~= 0
           B(n1,1)=B(n1,1)/B(n1,4);
           B(n1,2)=B(n1,2)/B(n1,4);
           B(n1,3)=B(n1,3)/B(n1,4);
        end

     end
    
toc;   

display(['  End TD calculation for E and B.']);

