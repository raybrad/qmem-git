function [Vs,ns,ps,Xs] = staticSolution(Vs,ns,ps,Xs,currdlink)

display('Start static solution under illum with QM');

global sigma epsilon_in epsilon_sd;
global mun mup;
global scl;
global nodes links contacts;
global Nnode;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes dirNodes;
global isBndNodes isDirNodes dcVolDirNodes;
global doping; 
global semiNodes;
global dirSemiNodes isDirSemiNodes;
global kx ky kz;
global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global sQMlinks;
global isqmnodes;
global bndLinks;
global metalNodes;
global ni;

global eqnNodes eqnSemiNodes;

global ehdifc;
global ehcrf;


%%%%%%% generate initial guess (equilibrium solution) %%%%%%%%%%%%%
eqnNodes = setdiff((1:Nnode)',[dirNodes;sQMnodes]);
eqnSemiNodes = setdiff(semiNodes,[dirSemiNodes;sQMnodes;dirNodes]);
eqnSemiNodesp  = setdiff(semiNodes,[sQMnodes]);
eqnSemiNodesp1 = unique([dirSemiNodes;dirNodes]);


%%% numbers of Poisson and DD equations
Nn = length(eqnNodes);
Nn_semi = length(eqnSemiNodes);
Nn_semip = length(eqnSemiNodesp);

%%%%%%%% start Newton's iteration %%%%%%%%
newtonTol = 1e-8;
updateTol = 1e-11;
maxNewtonIt = 15;
Eps = [epsilon_sd,epsilon_mt,epsilon_in];
Sgm = [0,sigma,0];
normRes = 1;
normUpdate = 1;
itNr = 0;
tic;

while itNr < maxNewtonIt && normUpdate > updateTol

    JF_v = zeros(Nnode,Nnode);  JF_n = zeros(Nnode,Nnode); 
    JF_p = zeros(Nnode,Nnode);  JF_x = zeros(Nnode,Nnode);  rhs_F = zeros(Nnode,1);

    JKn_v = zeros(Nnode,Nnode); JKn_n = zeros(Nnode,Nnode); 
    JKn_p = zeros(Nnode,Nnode); JKn_x = zeros(Nnode,Nnode); rhs_Kn = zeros(Nnode,1);  

    JKp_v = zeros(Nnode,Nnode); JKp_n = zeros(Nnode,Nnode);
    JKp_p = zeros(Nnode,Nnode); JKp_x = zeros(Nnode,Nnode); rhs_Kp = zeros(Nnode,1);

    JKx_v = zeros(Nnode,Nnode); JKx_n = zeros(Nnode,Nnode); 
    JKx_p = zeros(Nnode,Nnode); JKx_x = zeros(Nnode,Nnode); rhs_Kx = zeros(Nnode,1);

    %%%%%%%% build Jacobian and rhs of F %%%%%%%%%%%%%
    for n1 = eqnNodes.'
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        sign_n1 = sign(ajnd_n1-n1);

        if any(ajvolM_n1 == 1)  % node attached with semiconductor
            dV1 = deltaV(doping(n1));
            if all(ajvolM_n1 == 1) % bulk of semiconductor
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    n3 = links(lk,3);
                    if isSQMlinks(lk) == 0
                      ajvol_lk = linkVolumes{lk}(1,:);
                      ajvolS_lk = linkVolumes{lk}(2,:);
                      ajvolM_lk = volumeM(ajvol_lk);
                      coefV = sum(Eps(ajvolM_lk).*ajvolS_lk)/linkL(lk);
                      JF_v(n1,n1) = JF_v(n1,n1)+coefV;
                      JF_v(n1,n2) = JF_v(n1,n2)-coefV;
                      rhs_F(n1)   = rhs_F(n1) + epsilon_sd*linkS(lk)*(Vs(n1)-Vs(n2))/linkL(lk);
                    else
                      rhs_F(n1)   = rhs_F(n1) + 0;
                    end
                end
                rhs_F(n1) = rhs_F(n1) + (ns(n1)-ps(n1)-doping(n1))*nodeV(n1);
                JF_n(n1,n1) = nodeV(n1);
                JF_p(n1,n1) = -nodeV(n1);
           
            elseif any(ajvolM_n1 == 2)  % semi-metal interface or triple points
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    n3 = links(lk,3);
                    if isSQMlinks(lk) == 0
                       ajvol_lk = linkVolumes{lk}(1,:);
                       ajvolS_lk = linkVolumes{lk}(2,:);
                       ajvolM_lk = volumeM(ajvol_lk);
                       alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                       if isDirSemiNodes(n2)
                          dV2 = deltaV(doping(n2));
                          ttt = 0;
                       else
                          dV2 = 0;
                          ttt = 1;
                       end
                       for j = 1:length(ajvol_lk)
                           switch ajvolM_lk(j)
                                 case 1
                                      beta = Vs(n2)-Vs(n1);
                                      J = Jc(alpha_n,beta,'n',ns(n1),ns(n2))...
                                          +Jc(alpha_p,beta,'p',ps(n1),ps(n2));
                                      dJv = dJndV(alpha_n,beta,ns(n1),ns(n2))...
                                           +dJpdV(alpha_p,beta,ps(n1),ps(n2));
                                      dJn_nj = dJdcj(alpha_n,beta,'n');
                                      JF_n(n1,n2) = JF_n(n1,n2)+ ttt*ajvolS_lk(j)*dJn_nj;
                                      dJp_pj = dJdcj(alpha_p,beta,'p');
                                      JF_p(n1,n2) = JF_p(n1,n2)+ ttt*ajvolS_lk(j)*dJp_pj;
                                 case 2
                                      J = sigma*(Vs(n1)-Vs(n2)+dV1-dV2)/linkL(lk);
                                      dJv = sigma/linkL(lk);
                                 case 3
                                      J = 0;
                                      dJv = 0;
                                 otherwise
                                      error('undefined material');
                           end

                           JF_v(n1,n1) = JF_v(n1,n1)+ajvolS_lk(j)*dJv;
                           JF_v(n1,n2) = JF_v(n1,n2)-ajvolS_lk(j)*dJv;
                           rhs_F(n1)   = rhs_F(n1)+ajvolS_lk(j)*J;

                       end
                    else
                       rhs_F(n1)   = rhs_F(n1) + sign_n1(i)*linkS(lk)*currdlink(lk,n3+3);
                    end
                end                 
            elseif any(ajvolM_n1 == 3) % semi-insulator interface
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    n3 = links(lk,3);
                    if isSQMlinks(lk) == 0
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
                     JF_v(n1,n2) = JF_v(n1,n2)-coefV;
                     dV2_vec = [0,1,1/2];
                     rhs_F(n1) = rhs_F(n1)+sum(Eps(ajvolM_lk).*...
                        (-(Vs(n2)-Vs(n1)+dV2_vec(ajvolM_lk).*dV2)).*ajvolS_lk)/linkL(lk);
                    else
                     rhs_F(n1)   = rhs_F(n1) + 0;
                    end
                end
                semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
                JF_n(n1,n1) = semiV;
                JF_p(n1,n1) = -semiV;
                rhs_F(n1) = rhs_F(n1)+(ns(n1)-ps(n1)-doping(n1))*semiV;
            else
                error('Incorrect node assignment');
            end
        else                   % node attached no semiconductor
            for i = 1:length(ajlk_n1)
                n2 = ajnd_n1(i);
                lk = ajlk_n1(i);
                n3 = links(lk,3);
                if isSQMlinks(lk) == 0
                 if isDirSemiNodes(n2)
                    dV2 = deltaV(doping(n2));
                 else
                    dV2 = 0;
                 end
                 ajvol_lk = linkVolumes{lk}(1,:);
                 ajvolS_lk = linkVolumes{lk}(2,:);
                 ajvolM_lk = volumeM(ajvol_lk);

                 if any(ajvolM_n1 == 2)
                  tcoem = 0;
                 else
                  tcoem = 1;
                 end

                 coefV = sum((tcoem*Eps(ajvolM_lk)+Sgm(ajvolM_lk)).*ajvolS_lk)/linkL(lk);
                 JF_v(n1,n1) = JF_v(n1,n1)+coefV;
                 JF_v(n1,n2) = JF_v(n1,n2)-coefV;
                 dV2_vec = [0,1,1/2];
                 rhs_F(n1) = rhs_F(n1)+sum((tcoem*Eps(ajvolM_lk)+Sgm(ajvolM_lk)).*...
                        (-(Vs(n2)-Vs(n1)+dV2_vec(ajvolM_lk).*dV2)).*ajvolS_lk)/linkL(lk);
                elseif any(ajvolM_n1 == 2)
                  rhs_F(n1)   = rhs_F(n1) + sign_n1(i)*linkS(lk)*currdlink(lk,n3+3);
                else
                  rhs_F(n1)   = rhs_F(n1) + 0;
                end
            end
        end
        if any(ajvolM_n1 == 2) % scaling some equations with condutivity for matrix balance
            JF_v(n1,:) = JF_v(n1,:)/sigma;
            JF_n(n1,:) = JF_n(n1,:)/sigma;
            JF_p(n1,:) = JF_p(n1,:)/sigma;
            rhs_F(n1) = rhs_F(n1)/sigma;
        end
    end

    %%%%%%%%% build Jacobian and rhs of K %%%%%%%%%%%%%
    if ~isempty(eqnSemiNodes)
        chgngrns = chgngr(Vs,ns,ps,Xs);
        dchgngrdns = dchgngrd(Vs,ns,ps);
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
                n3 = links(lk,3);
                if isSQMlinks(lk) == 0
                 ajvol_lk = linkVolumes{lk}(1,:);
                 ajvolS_lk = linkVolumes{lk}(2,:);
                 ajvolM_lk = volumeM(ajvol_lk);
                 semiS = sum(ajvolS_lk(ajvolM_lk == 1));
                 alpha_n = mun/linkL(lk);  alpha_p = mup/linkL(lk);
                 beta = Vs(n2)-Vs(n1);
                 betax = (Xs(n2)-Xs(n1))/linkL(lk); 
 
                 Jx    =  ehdifc * betax;
                 dJx_x = ehdifc/linkL(lk);
                 rhs_Kx(n1) = rhs_Kx(n1)+Jx*semiS;
                 JKx_x(n1,n1) = JKx_x(n1,n1)-semiS*dJx_x;
                 JKx_x(n1,n2) = JKx_x(n1,n2)+semiS*dJx_x;
                 
                 dJn_v = dJndV(alpha_n,beta,ns(n1),ns(n2));
                 JKn_v(n1,n1) = JKn_v(n1,n1)+semiS*dJn_v;
                 JKn_v(n1,n2) = JKn_v(n1,n2)-semiS*dJn_v;
                 dJp_v = dJpdV(alpha_p,beta,ps(n1),ps(n2));
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
                
                 Jn = Jc(alpha_n,beta,'n',ns(n1),ns(n2));
                 Jp = Jc(alpha_p,beta,'p',ps(n1),ps(n2));
                 rhs_Kn(n1) = rhs_Kn(n1)+Jn*semiS;
                 rhs_Kp(n1) = rhs_Kp(n1)+Jp*semiS;
                else
                 rhs_Kn(n1)   = rhs_Kn(n1) + ehcrf     * sign_n1(i)*linkS(lk)*currdlink(lk,n3+3);
                 rhs_Kp(n1)   = rhs_Kp(n1) + (1-ehcrf) * sign_n1(i)*linkS(lk)*currdlink(lk,n3+3);
                 rhs_Kx(n1)   = rhs_Kx(n1) + 0;
                end
            end
            semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
            JKn_n(n1,n1) = JKn_n(n1,n1) + semiV*dchgngrdns(n1,1);
            JKn_p(n1,n1) = JKn_p(n1,n1) + semiV*dchgngrdns(n1,2);
            JKn_x(n1,n1) = JKn_x(n1,n1) + semiV*dchgngrdns(n1,3);
            JKp_p(n1,n1) = JKp_p(n1,n1) - semiV*dchgngrdns(n1,2);
            JKp_n(n1,n1) = JKp_n(n1,n1) - semiV*dchgngrdns(n1,1);
            JKp_x(n1,n1) = JKp_x(n1,n1) - semiV*dchgngrdns(n1,3);

            JKx_p(n1,n1) = JKx_p(n1,n1) - semiV*dchgngrdns(n1,2);
            JKx_n(n1,n1) = JKx_n(n1,n1) - semiV*dchgngrdns(n1,1);
            JKx_x(n1,n1) = JKx_x(n1,n1) + semiV*dchgngrdns(n1,4);

            rhs_Kn(n1) = rhs_Kn(n1) + semiV*chgngrns(n1,1);
            rhs_Kp(n1) = rhs_Kp(n1) - semiV*chgngrns(n1,1);
            rhs_Kx(n1) = rhs_Kx(n1) + semiV*chgngrns(n1,2);
        end
        for n1 = eqnSemiNodesp1.'
            ajlk_n1 = nodeLinks{n1}(1,:);
            ajnd_n1 = nodeLinks{n1}(2,:);
            ajvol_n1 = nodeVolumes{n1}(1,:);
            ajvolV_n1 = nodeVolumes{n1}(2,:);
            ajvolM_n1 = volumeM(ajvol_n1);
            sign_n1 = sign(ajnd_n1-n1);
            for i = 1:length(ajlk_n1)
                n2 = ajnd_n1(i);
                lk = ajlk_n1(i);
                n3 = links(lk,3);
                if isSQMlinks(lk) == 0
                  ajvol_lk = linkVolumes{lk}(1,:);
                  ajvolS_lk = linkVolumes{lk}(2,:);
                  ajvolM_lk = volumeM(ajvol_lk);
                  semiS = sum(ajvolS_lk(ajvolM_lk == 1));
                  betax = (Xs(n2)-Xs(n1))/linkL(lk);

                  Jx    =  ehdifc * betax;
                  dJx_x = ehdifc/linkL(lk);
                  rhs_Kx(n1) = rhs_Kx(n1)+Jx*semiS;
                  JKx_x(n1,n1) = JKx_x(n1,n1)-semiS*dJx_x;
                  JKx_x(n1,n2) = JKx_x(n1,n2)+semiS*dJx_x;
                else
                 rhs_Kx(n1)   = rhs_Kx(n1) + 0;
                end
            end
            semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
            JKx_p(n1,n1) = JKx_p(n1,n1) - semiV*dchgngrdns(n1,2);
            JKx_n(n1,n1) = JKx_n(n1,n1) - semiV*dchgngrdns(n1,1);
            JKx_x(n1,n1) = JKx_x(n1,n1) + semiV*dchgngrdns(n1,4);

            rhs_Kx(n1) = rhs_Kx(n1) + semiV*chgngrns(n1,2);
        end
    end

    JF_v = sparse(JF_v(eqnNodes,eqnNodes));
    JF_n = sparse(JF_n(eqnNodes,eqnSemiNodes));
    JF_p = sparse(JF_p(eqnNodes,eqnSemiNodes));
    JF_x = sparse(JF_x(eqnNodes,eqnSemiNodesp));

    JKn_v = sparse(JKn_v(eqnSemiNodes,eqnNodes));
    JKn_n = sparse(JKn_n(eqnSemiNodes,eqnSemiNodes));
    JKn_p = sparse(JKn_p(eqnSemiNodes,eqnSemiNodes));
    JKn_x = sparse(JKn_x(eqnSemiNodes,eqnSemiNodesp));

    JKp_v = sparse(JKp_v(eqnSemiNodes,eqnNodes));
    JKp_n = sparse(JKp_n(eqnSemiNodes,eqnSemiNodes));
    JKp_p = sparse(JKp_p(eqnSemiNodes,eqnSemiNodes));
    JKp_x = sparse(JKp_x(eqnSemiNodes,eqnSemiNodesp));

    JKx_v = sparse(JKx_v(eqnSemiNodesp,eqnNodes));
    JKx_n = sparse(JKx_n(eqnSemiNodesp,eqnSemiNodes));
    JKx_p = sparse(JKx_p(eqnSemiNodesp,eqnSemiNodes));
    JKx_x = sparse(JKx_x(eqnSemiNodesp,eqnSemiNodesp));

    rhs_F = sparse(rhs_F(eqnNodes));
    rhs_Kn = sparse(rhs_Kn(eqnSemiNodes));
    rhs_Kp = sparse(rhs_Kp(eqnSemiNodes));
    rhs_Kx = sparse(rhs_Kx(eqnSemiNodesp));

    Jacob = [JF_v,JF_n,JF_p,JF_x;JKn_v,JKn_n,JKn_p,JKn_x;JKp_v,JKp_n,JKp_p,JKp_x;JKx_v,JKx_n,JKx_p,JKx_x];
    rhs = -[rhs_F;rhs_Kn;rhs_Kp;rhs_Kx];
%    invJacob = Jacob^(-1);
%    dX = invJacob*rhs;
    dX = mylusovle(Jacob,rhs,1);
    dV = dX(1:Nn);
    dn = dX(Nn+1:Nn+Nn_semi);
    dp = dX(Nn+Nn_semi+1:Nn+2*Nn_semi);
    dxd= dX(Nn+2*Nn_semi+1:Nn+2*Nn_semi+Nn_semip);

    normRes = norm(rhs);
%    normUpdate = max(norm(dV)/norm(Vs(eqnNodes)));
    itNr = itNr+1;
    
    Vs(eqnNodes) = Vs(eqnNodes)+dV;
    ns(eqnSemiNodes) = ns(eqnSemiNodes)+dn;
    ps(eqnSemiNodes) = ps(eqnSemiNodes)+dp;    
    Xs(eqnSemiNodesp) = Xs(eqnSemiNodesp)+dxd;

%    Xs = excitonbc(Vs,ns,ps,Xs);

%    currb=currentS([2,round(ky/2),round(kz/2)],Vs,ns,ps,1);

%    currt=currentS([1,round(ky/2),round(kz/2)],Vs,ns,ps,1)-currentS([kx,round(ky/2),round(kz/2)],Vs,ns,ps,1)+currentS([round(kx/2),1,round(kz/2)],Vs,ns,ps,2)-currentS([round(kx/2),ky,round(kz/2)],Vs,ns,ps,2)+currentS([round(kx/2),round(ky/2),1],Vs,ns,ps,3)-currentS([round(kx/2),round(ky/2),kz],Vs,ns,ps,3);

%   currt=0;
	
%	normUpdate = max(norm(dV)/norm(Vs(eqnNodes)),norm(currb-curra)/norm(curra));

%    normUpdate = max([normRes,norm(dV)/norm(Vs(eqnNodes)),abs(currt)/abs(currb)]);
    normUpdate = max([normRes,norm(dV)/norm(Vs(eqnNodes)),norm(dn)/norm(ns(eqnSemiNodes)),norm(dp)/norm(ps(eqnSemiNodes)),norm(dxd)/norm(Xs(eqnSemiNodesp))]);

    display(['Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
	        '; ','normUpdate:',num2str(normUpdate),'; ','normdV:',num2str(norm(dV))]);
end  
% 
% Vs = Vs.*scl.Vt;
% ns = ns.*scl.ni;
% ps = ps.*scl.ni;

toc;   

display('End static solution under illum with QM');
