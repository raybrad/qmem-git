function [V,n,p,A,H] = tdrelaxpstep(V_p,n_p,p_p,A_p,H_p,dtV,dtn,dtp,dtH,dt)

display(['  Start TD simulation relaxation per step:']);

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
global doping; 
global semiNodes;
global dirSemiNodes isDirSemiNodes;
global kx ky kz;
global bndLinks;
global metalNodes;
global savefile;

%%%%%%% initial guess  %%%%%%%%%%%%%
V = V_p + dtV *dt;
n = n_p + dtn *dt;
p = p_p + dtp *dt;
H = H_p + dtH *dt;
A = A_p + H_p *dt + dtH * dt^2/2;


eqnNodes = setdiff((1:Nnode)',dirNodes);
eqnSemiNodes = setdiff(semiNodes,[dirSemiNodes;dirNodes]);
eqnLinks = setdiff((1:Nlink)',bndLinks);

%%%%%%%% start Newton's iteration %%%%%%%%
newtonTol = 1e-15;
updateTol = 1e-8;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 20;
maxLinIt = 200;
Eps = [epsilon_sd,epsilon_mt,epsilon_in];
Sgm = [0,sigma,0];

NeqnNodes = length(eqnNodes);
NeqnSemiNodes = length(eqnSemiNodes);
Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
normRes_pre = 1e10;
itNr = 0;
tic;
    
while itNr < maxNewtonIt && normUpdate > updateTol

     JF_v = zeros(Nnode,Nnode); JF_n = zeros(Nnode,Nnode);
     JF_p = zeros(Nnode,Nnode);
     JF_A = zeros(Nnode,Nlink); JF_H = zeros(Nnode,Nlink);
     JKn_v = zeros(Nnode,Nnode); JKn_n = zeros(Nnode,Nnode);
     JKn_p = zeros(Nnode,Nnode);
     JKn_A = zeros(Nnode,Nlink); JKn_H = zeros(Nnode,Nlink);
     JKp_v = zeros(Nnode,Nnode); JKp_n = zeros(Nnode,Nnode);
     JKp_p = zeros(Nnode,Nnode); JKp_A = zeros(Nnode,Nlink); JKp_H = zeros(Nnode,Nlink);
     JG_v = zeros(Nlink,Nnode); JG_n = zeros(Nlink,Nnode);
     JG_p = zeros(Nlink,Nnode);
     JG_A = zeros(Nlink,Nlink); JG_H = zeros(Nlink,Nlink);
     
     JT_v = zeros(Nlink,Nnode); JT_n = zeros(Nlink,Nnode);
     JT_p = zeros(Nlink,Nnode);
     JT_A = zeros(Nlink,Nlink); JT_H = zeros(Nlink,Nlink);
          
     rhs_F = zeros(Nnode,1); rhs_Kn = zeros(Nnode,1);
     rhs_Kp = zeros(Nnode,1);
     rhs_G = zeros(Nlink,1);
     rhs_T = zeros(Nlink,1);
     
     for n1 = eqnNodes.'
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        sign_n1 = sign(ajnd_n1-n1);

        if isDirSemiNodes(n1), dV1 = deltaV(doping(n1)); else dV1 = 0; end

        if any(ajvolM_n1 == 2) % metal/semi interfaces or triple points
            for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    ajvol_lk = linkVolumes{lk}(1,:);
                    ajvolS_lk = linkVolumes{lk}(2,:);
                    ajvolM_lk = volumeM(ajvol_lk);
 
                    if isDirSemiNodes(n2)
                      dV2 = deltaV(doping(n2)); 
                      ttt = 0;
                    else 
                      dV2 = 0;
                      ttt = 1; 
                    end

                    dtE1 = -(dtV(n2)-dtV(n1))/linkL(lk)-sign_n1(i)*dtH(lk);
                    for j = 1:length(ajvol_lk)
                        switch ajvolM_lk(j)
                            case 1
                                alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                                beta = V(n2)-V(n1) + sign_n1(i)*H(lk)*linkL(lk);
                                J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                                    +Jc(alpha_p,beta,'p',p(n1),p(n2));
                                Jd = epsilon_sd * dtE1;
                                J = J_diff +Jd;
                                dJv = dJndV(alpha_n,beta,n(n1),n(n2))...
                                    +dJpdV(alpha_p,beta,p(n1),p(n2)) + epsilon_sd/(dt*linkL(lk));
                                dJH = -sign_n1(i)*dJv*linkL(lk);
                                dJn_nj = dJdcj(alpha_n,beta,'n');
                                JF_n(n1,n2) = JF_n(n1,n2)+ttt*ajvolS_lk(j)*dJn_nj;
                                dJp_pj = dJdcj(alpha_p,beta,'p');
                                JF_p(n1,n2) = JF_p(n1,n2)+ttt*ajvolS_lk(j)*dJp_pj;
                            case 2
                                J =  sigma * ((V(n1)-V(n2)+dV1-dV2)/linkL(lk) - sign_n1(i) * H(lk)) + epsilon_mt*dtE1;
                                dJv = sigma/linkL(lk) + epsilon_mt/(dt*linkL(lk));
                                dJH = -sign_n1(i)* (sigma+epsilon_mt/dt);
                            case 3
                                J   = epsilon_in*dtE1;
                                dJv = epsilon_in/(dt*linkL(lk));
                                dJH = -sign_n1(i)*epsilon_in/dt;
                            otherwise
                                error('undefined material');
                        end
                        JF_v(n1,n1) = JF_v(n1,n1)+ajvolS_lk(j)*dJv;
                        JF_v(n1,n2) = JF_v(n1,n2)-ajvolS_lk(j)*dJv;
                        JF_H(n1,lk) = JF_H(n1,lk)+ajvolS_lk(j)*dJH;
                        rhs_F(n1) = rhs_F(n1)+ajvolS_lk(j)*J;
                    end
            end
        else
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    if isDirSemiNodes(n2)
                        dV2 = deltaV(doping(n2));
                    else
                        dV2 = 0;
                    end
                    ajvol_lk = linkVolumes{lk}(1,:);
                    ajvolS_lk = linkVolumes{lk}(2,:);
                    ajvolM_lk = volumeM(ajvol_lk);
                    coefV = sum(Eps(ajvolM_lk).*ajvolS_lk)/linkL(lk);
                    JF_v(n1,n1) = JF_v(n1,n1)+coefV;
                    JF_v(n1,n2) = -coefV;
                    dV2_vec = [0,1,1/2];
                    coefH = -coefV*linkL(lk);
                    JF_H(n1,lk) = sign_n1(i)*coefH;
                    rhs_F(n1) = rhs_F(n1)+sum(Eps(ajvolM_lk).*(-(V(n2)-V(n1)+dV2_vec(ajvolM_lk).*dV2)/linkL(lk)-sign_n1(i)* H(lk)).*ajvolS_lk);
                end
                semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
                JF_n(n1,n1) = semiV;
                JF_p(n1,n1) = -semiV;
                rhs_F(n1) = rhs_F(n1)+(n(n1)-p(n1) - doping(n1))*semiV;
        end
        if any(ajvolM_n1 == 2) % scaling some equations with condutivity for matrix balance
            JF_v(n1,:) = JF_v(n1,:)/sigma;
            JF_n(n1,:) = JF_n(n1,:)/sigma;
            JF_p(n1,:) = JF_p(n1,:)/sigma;
            JF_A(n1,:) = JF_A(n1,:)/sigma;
            JF_H(n1,:) = JF_H(n1,:)/sigma;
            rhs_F(n1) = rhs_F(n1)/sigma;
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
                E1 = -(V(n2)-V(n1))/linkL(lk)-sign_n1(i)*H(lk);
                beta = - linkL(lk) * E1;

                dJn_v = dJndV(alpha_n,beta,n(n1),n(n2));
                JKn_v(n1,n1) = JKn_v(n1,n1)+semiS*dJn_v;
                JKn_v(n1,n2) = JKn_v(n1,n2)-semiS*dJn_v;
                dJp_v = dJpdV(alpha_p,beta,p(n1),p(n2));
                JKp_v(n1,n1) = JKp_v(n1,n1)+semiS*dJp_v;
                JKp_v(n1,n2) = JKp_v(n1,n2)-semiS*dJp_v;

                dJn_ni = dJdci(alpha_n,beta,'n');
                dJn_nj = dJdcj(alpha_n,beta,'n');
                JKn_n(n1,n1) = JKn_n(n1,n1)+semiS*dJn_ni;
                JKn_n(n1,n2) = JKn_n(n1,n2)+semiS*dJn_nj;
                dJp_pi = dJdci(alpha_p,beta,'p');
                dJp_pj = dJdcj(alpha_p,beta,'p');
                JKp_p(n1,n1) = JKp_p(n1,n1)+semiS*dJp_pi;
                JKp_p(n1,n2) = JKp_p(n1,n2)+semiS*dJp_pj;

                JKn_H(n1,lk) = JKn_H(n1,lk)-sign_n1(i)*semiS*dJn_v*linkL(lk);
                JKp_H(n1,lk) = JKp_H(n1,lk)-sign_n1(i)*semiS*dJp_v*linkL(lk);

                Jn = Jc(alpha_n,beta,'n',n(n1),n(n2));
                Jp = Jc(alpha_p,beta,'p',p(n1),p(n2));
                rhs_Kn(n1) = rhs_Kn(n1)+Jn*semiS;
                rhs_Kp(n1) = rhs_Kp(n1)+Jp*semiS;
            end
            semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
            JKn_n(n1,n1) = JKn_n(n1,n1)-semiV/dt;
            JKp_p(n1,n1) = JKp_p(n1,n1)+semiV/dt;
            rhs_Kn(n1) = rhs_Kn(n1)-semiV*dtn(n1);
            rhs_Kp(n1) = rhs_Kp(n1)+semiV*dtp(n1);
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
        GD = zeros(1,Nlink);

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
            rhs_G(l1) = rhs_G(l1)+CC(ajlk)*A(ajlk);
        end
        rhs_G(l1) = rhs_G(l1)+CC(l1)*A(l1);

        %%%%%%%% gradient divergence of A and gradient of V %%%%%%%%
        if ~isBndNodes(n1)
            ajlk_n1 = nodeLinks{n1}(1,:);
            ajnd_n1 = nodeLinks{n1}(2,:);
            sign_n1 = sign(ajnd_n1-n1);
            coef_n1 = linkS(ajlk_n1)./nodeV(n1)*linkS(l1)/linkL(l1);
            GD(ajlk_n1) = GD(ajlk_n1)-sign_n1.*coef_n1';
            ajvol_n1 = nodeVolumes{n1}(1,:);
            ajvolV_n1 = nodeVolumes{n1}(2,:);
            ajvolM_n1 = volumeM(ajvol_n1);
            coefV_n1 = scl.K*(sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1))*linkS(l1)/linkL(l1);

            JG_v(l1,n1) = JG_v(l1,n1)+coefV_n1/dt;

            rhs_G(l1) = rhs_G(l1)-(-sign_n1.*coef_n1')*A(ajlk_n1);
            rhs_G(l1) = rhs_G(l1)+coefV_n1*dtV(n1);
        end

        if ~isBndNodes(n2)
            ajlk_n2 = nodeLinks{n2}(1,:);
            ajnd_n2 = nodeLinks{n2}(2,:);
            sign_n2 = sign(ajnd_n2-n2);
            coef_n2 = linkS(ajlk_n2)./nodeV(n2)*linkS(l1)/linkL(l1);
            GD(ajlk_n2) = GD(ajlk_n2)+sign_n2.*coef_n2';
            ajvol_n2 = nodeVolumes{n2}(1,:);
            ajvolV_n2 = nodeVolumes{n2}(2,:);
            ajvolM_n2 = volumeM(ajvol_n2);
            coefV_n2 = scl.K*(sum(Eps(ajvolM_n2).*ajvolV_n2)/nodeV(n2))*linkS(l1)/linkL(l1);

            JG_v(l1,n2) = JG_v(l1,n2)-coefV_n2/dt;

            rhs_G(l1) = rhs_G(l1)-sign_n2.*coef_n2'*A(ajlk_n2);
            rhs_G(l1) = rhs_G(l1)-coefV_n2*dtV(n2);
        end

        JG_A(l1,:) = CC-GD;

        %%%%%% source current %%%%%%%%%%%
        E1 = -(V(n2)-V(n1))/linkL(l1)-H(l1);
        dtE1 = -(dtV(n2)-dtV(n1))/linkL(l1) - dtH(l1);
        if isDirSemiNodes(n1), dV1 = deltaV(doping(n1)); else dV1 = 0; end
        if isDirSemiNodes(n2), dV2 = deltaV(doping(n2)); else dV2 = 0; end
        for i = 1:length(ajvol_l1)
            switch ajvolM_l1(i)
                case 1
                    alpha_n = mun/linkL(l1); alpha_p = mup/linkL(l1);
                    beta = V(n2)-V(n1)+H(l1)*linkL(l1);
                    J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                        +Jc(alpha_p,beta,'p',p(n1),p(n2));
                    J = J_diff + epsilon_sd*dtE1;
                    dJ_v = dJndV(alpha_n,beta,n(n1),n(n2))...
                          +dJpdV(alpha_p,beta,p(n1),p(n2))+epsilon_sd/(dt*linkL(l1));
                    dJn_ni = dJdci(alpha_n,beta,'n');
                    dJn_nj = dJdcj(alpha_n,beta,'n');
                    dJp_pi = dJdci(alpha_p,beta,'p');
                    dJp_pj = dJdcj(alpha_p,beta,'p');
                    dJ_H = - dJ_v*linkL(l1); 
                    JG_n(l1,n1) = JG_n(l1,n1)-scl.K*ajvolS_l1(i)*dJn_ni;
                    JG_n(l1,n2) = JG_n(l1,n2)-scl.K*ajvolS_l1(i)*dJn_nj;
                    JG_p(l1,n1) = JG_p(l1,n1)-scl.K*ajvolS_l1(i)*dJp_pi;
                    JG_p(l1,n2) = JG_p(l1,n2)-scl.K*ajvolS_l1(i)*dJp_pj;
                case 2
                    dJ_v = sigma/linkL(l1)+epsilon_mt/(dt*linkL(l1));
                    dJ_H = -(sigma+epsilon_mt/dt);
                    J = sigma*((V(n1)-V(n2)+dV1-dV2)/linkL(l1) - H(l1)) + epsilon_mt*dtE1;
                case 3
                    dJ_v = epsilon_in/(dt*linkL(l1));
                    dJ_H = -epsilon_in/dt;
                    J = epsilon_in*dtE1;
                otherwise
                    error('undefined material');
            end
            JG_v(l1,n1) = JG_v(l1,n1)-scl.K*ajvolS_l1(i)*dJ_v;
            JG_v(l1,n2) = JG_v(l1,n2)+scl.K*ajvolS_l1(i)*dJ_v;
            JG_H(l1,l1) = JG_H(l1,l1)-scl.K*ajvolS_l1(i)*dJ_H;
            rhs_G(l1) = rhs_G(l1)-scl.K*ajvolS_l1(i)*J;
        end
        
        rhs_T(l1)  = H(l1)*dt - dtH(l1) * dt^2/2 + A_p(l1) - A(l1);
        JT_A(l1,l1)= -1;
        JT_H(l1,l1)= dt/2;
     end
    
    JF_v = JF_v(eqnNodes,eqnNodes);
    JF_n = JF_n(eqnNodes,eqnSemiNodes);
    JF_p = JF_p(eqnNodes,eqnSemiNodes);
    JF_A  = JF_A(eqnNodes,eqnLinks);
    JF_H  = JF_H(eqnNodes,eqnLinks);
    
    JKn_v = JKn_v(eqnSemiNodes,eqnNodes);
    JKn_n = JKn_n(eqnSemiNodes,eqnSemiNodes);
    JKn_p = JKn_p(eqnSemiNodes,eqnSemiNodes);
    JKn_A  = JKn_A(eqnSemiNodes,eqnLinks);
    JKn_H  = JKn_H(eqnSemiNodes,eqnLinks);
    
    JKp_v = JKp_v(eqnSemiNodes,eqnNodes);
    JKp_n = JKp_n(eqnSemiNodes,eqnSemiNodes);
    JKp_p = JKp_p(eqnSemiNodes,eqnSemiNodes);
    JKp_A  = JKp_A(eqnSemiNodes,eqnLinks);
    JKp_H  = JKp_H(eqnSemiNodes,eqnLinks);
    
    JG_v = JG_v(eqnLinks,eqnNodes);
    JG_n = JG_n(eqnLinks,eqnSemiNodes);
    JG_p = JG_p(eqnLinks,eqnSemiNodes);
    JG_A  = JG_A(eqnLinks,eqnLinks);
    JG_H  = JG_H(eqnLinks,eqnLinks);
    
    JT_v = JT_v(eqnLinks,eqnNodes);
    JT_n = JT_n(eqnLinks,eqnSemiNodes);
    JT_p = JT_p(eqnLinks,eqnSemiNodes);
    JT_A  = JT_A(eqnLinks,eqnLinks);
    JT_H  = JT_H(eqnLinks,eqnLinks);
    
    rhs_F = rhs_F(eqnNodes);
    rhs_Kn = rhs_Kn(eqnSemiNodes);
    rhs_Kp = rhs_Kp(eqnSemiNodes);
    rhs_G = rhs_G(eqnLinks);
    rhs_T = rhs_T(eqnLinks);

    JF = [JF_v,JF_n,JF_p,JF_A,JF_H];
    JKn = [JKn_v,JKn_n,JKn_p,JKn_A,JKn_H];
    JKp = [JKp_v,JKp_n,JKp_p,JKp_A,JKp_H];
    JG = [JG_v,JG_n,JG_p,JG_A,JG_H];
    JT = [JT_v,JT_n,JT_p,JT_A,JT_H];

    clear JF_v JF_n JF_p JF_A JF_H;
    clear JKn_v JKn_n JKn_p JKn_A JKn_H;
    clear JKp_v JKp_n JKp_p JKp_A JKp_H;
    clear JG_v JG_n JG_p JG_A JG_H;
      
    Jacob = [JF;JKn;JKp;JG;JT];
    rhs = -[rhs_F;rhs_Kn;rhs_Kp;rhs_G;rhs_T];
    
    clear JF JKn JKp JG JT;
    clear rhs_F rhs_Kn rhs_Kp rhs_G rhs_T;
%     Jacob1 = sparse(Jacob);
%     rhs1 = sparse(rhs);
%     [L,U] = luinc(Jacob1,droptol);
%     [dX,flag,relres,iter] = gmres(Jacob1,rhs1,[],linSolveTol,maxLinIt,L,U);
%     if flag ~= 0, error('No convergence for linear solver'); end
%     gmresItNr(l) = iter(2);
%    invJacob = Jacob^(-1);
%    dX = invJacob*rhs;
%    dX = Jacob\rhs;
    dX = mylusovle(Jacob,rhs,2);
    dV = dX(1:NeqnNodes);
    dn = dX(NeqnNodes+1:NeqnNodes+NeqnSemiNodes);
    dp = dX(NeqnNodes+NeqnSemiNodes+1:NeqnNodes+2*NeqnSemiNodes);
    dA = dX(NeqnNodes+2*NeqnSemiNodes+1:NeqnNodes+2*NeqnSemiNodes+Nl);
    dH = dX(NeqnNodes+2*NeqnSemiNodes+Nl+1:NeqnNodes+2*NeqnSemiNodes+2*Nl);
    
    normRes = norm(rhs);
    normRes_pre = normRes;

    itNr = itNr+1;
    
    V(eqnNodes) = V(eqnNodes)+dV;
    n(eqnSemiNodes) = n(eqnSemiNodes)+dn;
    p(eqnSemiNodes) = p(eqnSemiNodes)+dp;
    A(eqnLinks) = A(eqnLinks)+dA;
    H(eqnLinks) = H(eqnLinks)+dH;

    dtV  = (V - V_p)/dt;
    dtH  = (H - H_p)/dt;
    dtn  = (n - n_p)/dt;
    dtp  = (p - p_p)/dt;

    currt=currentStd([1,round(ky/2),round(kz/2)],V,n,p,H,dtV,dtH,1)-currentStd([kx,round(ky/2),round(kz/2)],V,n,p,H,dtV,dtH,1)+currentStd([round(kx/2),1,round(kz/2)],V,n,p,H,dtV,dtH,2)-currentStd([round(kx/2),ky,round(kz/2)],V,n,p,H,dtV,dtH,2)+currentStd([round(kx/2),round(ky/2),1],V,n,p,H,dtV,dtH,3)-currentStd([round(kx/2),round(ky/2),kz],V,n,p,H,dtV,dtH,3);

    currb = currentStd([2,round(ky/2),round(kz/2)],V,n,p,H,dtV,dtH,1);

    normUpdate = max([normRes,norm(dV)/norm(V(eqnNodes)),norm(dA)/norm(A(eqnLinks)),norm(dH)/norm(H(eqnLinks)),abs(currt)/abs(currb)]);

%    normUpdate = max([norm(dV)/norm(V(eqnNodes)),norm(dA)/norm(A(eqnLinks)),norm(dB)/norm(B(eqnLinks))]);


    display(['  Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
        '; ','normUpdate:',num2str(normUpdate),'; ','normdV:',num2str(norm(dV))]);

    clear dX dV dn dp dA dH;
    
end    
toc;   

display(['  End TD simulation relaxation per step.']);

