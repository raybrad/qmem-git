function [dtV,dtn,dtp,dtH] = tdcalupdates(V,n,p,A,H,dtVp,dtHp)

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
global isBndNodes isDirNodes dcVolDirNodes acVolDirNodes;
global doping; % doping profile
global semiNodes;
global dirSemiNodes isDirSemiNodes;
global kx ky kz;
global bndLinks;
global metalNodes;

%%%%%%% initial guess  %%%%%%%%%%%%%

eqnNodes  = setdiff((1:Nnode)',dirNodes);
eqnSemiNodes = setdiff(semiNodes,[dirSemiNodes;dirNodes]);
eqnLinks = setdiff((1:Nlink)',bndLinks);

%%%%%%%% start Newton's iteration %%%%%%%%
Eps     = [epsilon_sd,epsilon_mt,epsilon_in];
ntp  = 1;
tic;

dtV  = dtVp; 
dtn  = zeros(Nnode,1); 
dtp  = zeros(Nnode,1); 
dtH  = dtHp;
    
%%%%%%%% build Jacobian and rhs of F (Gauss's law) %%%%%%%%%%%%%
for n1 = eqnNodes.'
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        sign_n1 = sign(ajnd_n1-n1);
        if ~isBndNodes(n1)
          coef_n1 = linkS(ajlk_n1)./nodeV(n1);
          dtV(n1)  = -sign_n1.*coef_n1'*A(ajlk_n1);
          epsnode = sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1);
          dtV(n1) = dtV(n1) / (scl.K * epsnode);
        else
          eff = 0;
          if isDirSemiNodes(n1), dV1 = deltaV(doping(n1)); else dV1 = 0; end
          for i = 1:length(ajlk_n1)
            n2 = ajnd_n1(i);
            lk = ajlk_n1(i);
            ajvol_lk = linkVolumes{lk}(1,:);
            ajvolS_lk = linkVolumes{lk}(2,:);
            ajvolM_lk = volumeM(ajvol_lk);
            if isDirSemiNodes(n2), dV2 = deltaV(doping(n2)); else dV2 = 0; end

            dtE1 = -dtVp(n2)/linkL(lk)-sign_n1(i)*dtHp(lk);
            for j = 1:length(ajvol_lk)
               switch ajvolM_lk(j)
                 case 1
                  alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                  beta = V(n2)-V(n1) + sign_n1(i)*H(lk)*linkL(lk);
                  J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                          +Jc(alpha_p,beta,'p',p(n1),p(n2));
                  Jd = epsilon_sd * dtE1;
                  J = J_diff +Jd;
                 case 2
                  J =  sigma * ((V(n1)-V(n2)+dV1-dV2)/linkL(lk) - sign_n1(i) * H(lk)) + epsilon_mt*dtE1;
                 case 3
                  J   = epsilon_in*dtE1;
                 otherwise
                  error('undefined material');
               end
               dtV(n1) = dtV(n1)+ajvolS_lk(j)*J;
            end
            eff=eff+sum(Eps(ajvolM_lk).*ajvolS_lk)/linkL(lk);
          end
          if eff~=0;
           dtV(n1) = -dtV(n1)/eff;
          end
        end
end
     
    %%%%%%%%% build Jacobian and rhs of K (current continuity of semi) %%%%%%%%%%%%%
if ~isempty(eqnSemiNodes)
        for n1 = eqnSemiNodes.'
            ajlk_n1 = nodeLinks{n1}(1,:);
            ajnd_n1 = nodeLinks{n1}(2,:);
            ajvol_n1 = nodeVolumes{n1}(1,:);
            ajvolV_n1 = nodeVolumes{n1}(2,:);
            ajvolM_n1 = volumeM(ajvol_n1);
            sign_n1 = sign(ajnd_n1-n1);
            for i = 1:length(ajlk_n1)
                n2 = ajnd_n1(i);
                lk = ajlk_n1(i);
                ajvol_lk = linkVolumes{lk}(1,:);
                ajvolS_lk = linkVolumes{lk}(2,:);
                ajvolM_lk = volumeM(ajvol_lk);
                semiS = sum(ajvolS_lk(ajvolM_lk == 1));
                alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                beta = V(n2)-V(n1) + sign_n1(i)*H(lk)*linkL(lk);
                Jn = Jc(alpha_n,beta,'n',n(n1),n(n2));
                Jp = Jc(alpha_p,beta,'p',p(n1),p(n2));
                dtn(n1) = dtn(n1)+Jn*semiS;
                dtp(n1) = dtp(n1)-Jp*semiS;
            end
            semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
            if semiV ~=0
              dtn(n1) = dtn(n1)/semiV;
              dtp(n1) = dtp(n1)/semiV;
            end
        end
end
    
    %%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
for l1 = eqnLinks'
        n1 = links(l1,1);
        n2 = links(l1,2);
        ajvol_l1 = linkVolumes{l1}(1,:);
        ajvolS_l1 = linkVolumes{l1}(2,:);
        ajvolM_l1 = volumeM(ajvol_l1);
        CC = zeros(1,Nlink);

        epslinkf = scl.K * sum(Eps(ajvolM_l1).*ajvolS_l1);
        
        %%%%%%%%%% curl curl of A %%%%%%%%%%
        ajsf_l1 = linkSurfs{l1};
        for i = 1:size(ajsf_l1,2)
            ajlk = ajsf_l1(2:4,i);
            S = linkL(l1)*linkL(ajlk(1));
            dL = dlinkL(ajsf_l1(1,i));
            ajlkL = linkL(ajlk)';
            if n2 == links(ajlk(1),1)
                ajlkSign = [1,-1,-1];
            elseif n2 == links(ajlk(1),2)
                ajlkSign = [-1,-1,1];
            else
                error('Incorrect link arrangement');
            end
            CC(ajlk) = dL/S*(ajlkSign.*ajlkL);
            CC(l1) = CC(l1)+dL/S*linkL(l1);
            dtH(l1) = dtH(l1)+CC(ajlk)*A(ajlk);
        end

        dtH(l1) = dtH(l1) + CC(l1)*A(l1);
        
        %%%%%% source current %%%%%%%%%%%
        E1 = -(V(n2)-V(n1))/linkL(l1)-H(l1);    
        dtE1 = -(dtV(n2)-dtV(n1))/linkL(l1);
        if isDirSemiNodes(n1), dV1 = deltaV(doping(n1)); else dV1 = 0; end
        if isDirSemiNodes(n2), dV2 = deltaV(doping(n2)); else dV2 = 0; end
        for i = 1:length(ajvol_l1)
            switch ajvolM_l1(i)
                case 1
                    alpha_n = mun/linkL(l1); alpha_p = mup/linkL(l1);
                    beta = V(n2)-V(n1) + H(l1)*linkL(l1);
                    J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                        +Jc(alpha_p,beta,'p',p(n1),p(n2));
                    J = J_diff + epsilon_sd*dtE1;
                case 2
                    J = sigma*((V(n1)-V(n2)+dV1-dV2)/linkL(l1) - H(l1))+epsilon_mt*dtE1;
                case 3
                    J = epsilon_in*dtE1;
                otherwise
                    error('undefined material');
            end
            dtH(l1) = dtH(l1)-scl.K*ajvolS_l1(i)*J;
        end
     
        dtH(l1) = - dtH(l1)/epslinkf;
end
