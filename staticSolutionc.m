function [Vs,ns,ps] = staticSolution(Vs,ns,ps,currdlink)

%%% static solution of Poisson and drift-diffusion (DD) equations 
%%% use Newton's method
%%% only electron density n is used (p is all zero at this moment)

display('Start static solution cdbc');

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


%%% numbers of Poisson and DD equations
%dirSemiNodes = [dirSemiNodes;QMnodes];
eqnNodes = setdiff((1:Nnode)',[dirNodes;sQMnodes]);
eqnSemiNodes = setdiff(semiNodes,[dirSemiNodes;sQMnodes;dirNodes]);
Nn = length(eqnNodes);
Nn_semi = length(eqnSemiNodes);

%%%%%%%% start Newton's iteration %%%%%%%%
newtonTol = 1e-8;
updateTol = 1e-10;
maxNewtonIt = 20;
Eps = [epsilon_sd,epsilon_mt,epsilon_in];
Sgm = [0,sigma,0];
normRes = 1;
normUpdate = 1;
itNr = 0;
tic;

curra = ones(kx,1);

while itNr < maxNewtonIt && normUpdate > updateTol
    %tic;
    JF_v = zeros(Nnode,Nnode); JF_n = zeros(Nnode,Nnode); 
    JKn_v = zeros(Nnode,Nnode); JKn_n = zeros(Nnode,Nnode); 
    rhs_F = zeros(Nnode,1); rhs_Kn = zeros(Nnode,1); 
%    JF_p = []; JKn_p = []; JKp_v = []; JKp_n = [];JKp_p = []; rhs_Kp = [];
    JF_p = zeros(Nnode,Nnode); JKn_p = zeros(Nnode,Nnode); JKp_v = zeros(Nnode,Nnode);
    JKp_n = zeros(Nnode,Nnode);JKp_p = zeros(Nnode,Nnode); rhs_Kp = zeros(Nnode,1);
    currb=[];

    %%%%%%%% build Jacobian and rhs of F %%%%%%%%%%%%%
    for n1 = eqnNodes.'
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        if any(ajvolM_n1 == 1)  % node attached with semiconductor
            dV1 = deltaV(doping(n1));
            if all(ajvolM_n1 == 1) % bulk of semiconductor
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
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
                    elseif isSQMlinks(lk) == 1
                       rhs_F(n1)   = rhs_F(n1) + linkS(lk)*currdlink(lk,4);
                    elseif isSQMlinks(lk) == 2
                       rhs_F(n1)   = rhs_F(n1) - linkS(lk)*currdlink(lk,4);
                    elseif isSQMlinks(lk) == 3  
                       rhs_F(n1)   = rhs_F(n1) - linkS(lk)*currdlink(lk,6);
                    elseif isSQMlinks(lk) == 4
                       rhs_F(n1)   = rhs_F(n1) + 0;
                    elseif isSQMlinks(lk) > 6
                       rhs_F(n1)   = rhs_F(n1) + 0;
                    end
                end                 
            elseif any(ajvolM_n1 == 3) % semi-insulator interface
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
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
                elseif isSQMlinks(lk) == 1 && any(ajvolM_n1 == 2)
                  rhs_F(n1)   = rhs_F(n1) + linkS(lk)*currdlink(lk,4);
                elseif isSQMlinks(lk) == 2 && any(ajvolM_n1 == 2)
                  rhs_F(n1)   = rhs_F(n1) - linkS(lk)*currdlink(lk,4);
                elseif isSQMlinks(lk) == 3 && any(ajvolM_n1 == 2)
                  rhs_F(n1)   = rhs_F(n1) - linkS(lk)*currdlink(lk,6);
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
        for n1 = eqnSemiNodes.'
            ajlk_n1 = nodeLinks{n1}(1,:);
            ajnd_n1 = nodeLinks{n1}(2,:);
%             ajvol_n1 = nodeVolumes{n1}(1,:);
            for i = 1:length(ajlk_n1)
                n2 = ajnd_n1(i);
                lk = ajlk_n1(i);
                if isSQMlinks(lk) == 0
                 ajvol_lk = linkVolumes{lk}(1,:);
                 ajvolS_lk = linkVolumes{lk}(2,:);
                 ajvolM_lk = volumeM(ajvol_lk);
                 semiS = sum(ajvolS_lk(ajvolM_lk == 1));
                 alpha_n = mun/linkL(lk);  alpha_p = mup/linkL(lk);

                 if semiS ~= 0
                  beta = Vs(n2)-Vs(n1);
                
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
                 end
                elseif isSQMlinks(lk) == 1
                 rhs_Kn(n1) = rhs_Kn(n1) + ehcrf     * linkS(lk)*currdlink(lk,4);
                 rhs_Kp(n1) = rhs_Kp(n1) + (1-ehcrf) * linkS(lk)*currdlink(lk,4);
                elseif isSQMlinks(lk) == 2
                 rhs_Kn(n1) = rhs_Kn(n1) - ehcrf     * linkS(lk)*currdlink(lk,4);
                 rhs_Kp(n1) = rhs_Kp(n1) - (1-ehcrf) * linkS(lk)*currdlink(lk,4);
                elseif isSQMlinks(lk) == 3
                 rhs_Kn(n1) = rhs_Kn(n1) - ehcrf     * linkS(lk)*currdlink(lk,6);
                 rhs_Kp(n1) = rhs_Kp(n1) - (1-ehcrf) * linkS(lk)*currdlink(lk,6);
                else
                 rhs_Kn(n1) = rhs_Kn(n1) + 0;
                 rhs_Kp(n1) = rhs_Kp(n1) + 0;
                end
            end
        end
    end

    JF_v = sparse(JF_v(eqnNodes,eqnNodes));
    JF_n = sparse(JF_n(eqnNodes,eqnSemiNodes));
    JF_p = sparse(JF_p(eqnNodes,eqnSemiNodes));
    JKn_v = sparse(JKn_v(eqnSemiNodes,eqnNodes));
    JKn_n = sparse(JKn_n(eqnSemiNodes,eqnSemiNodes));
    JKn_p = sparse(JKn_p(eqnSemiNodes,eqnSemiNodes));
    JKp_v = sparse(JKp_v(eqnSemiNodes,eqnNodes));
    JKp_n = sparse(JKp_n(eqnSemiNodes,eqnSemiNodes));
    JKp_p = sparse(JKp_p(eqnSemiNodes,eqnSemiNodes));
    rhs_F = sparse(rhs_F(eqnNodes));
    rhs_Kn = sparse(rhs_Kn(eqnSemiNodes));
    rhs_Kp = sparse(rhs_Kp(eqnSemiNodes));
    Jacob = [JF_v,JF_n,JF_p;JKn_v,JKn_n,JKn_p;JKp_v,JKp_n,JKp_p];
    rhs = -[rhs_F;rhs_Kn;rhs_Kp];
%    invJacob = Jacob^(-1);
%    dX = invJacob*rhs;
    dX = mylusovle(Jacob,rhs,1);
    dV = dX(1:Nn);
    dn = dX(Nn+1:Nn+Nn_semi);
    dp = dX(Nn+Nn_semi+1:Nn+2*Nn_semi);
    normRes = norm(rhs);
%    normUpdate = max(norm(dV)/norm(Vs(eqnNodes)));
    itNr = itNr+1;
    
    Vs(eqnNodes) = Vs(eqnNodes)+dV;
    ns(eqnSemiNodes) = ns(eqnSemiNodes)+dn;
    ps(eqnSemiNodes) = ps(eqnSemiNodes)+dp;    

    for i=1:kx
        currb=[currb;currentS([i,round(ky/2),round(kz/2)],Vs,ns,ps,1)];
    end

%	currt=currentS([1,round(ky/2),round(kz/2)],Vs,ns,ps,1)-currentS([kx,round(ky/2),round(kz/2)],Vs,ns,ps,1)+currentS([round(kx/2),1,round(kz/2)],Vs,ns,ps,2)-currentS([round(kx/2),ky,round(kz/2)],Vs,ns,ps,2)+currentS([round(kx/2),round(ky/2),1],Vs,ns,ps,3)-currentS([round(kx/2),round(ky/2),kz],Vs,ns,ps,3);


    currt = 0;

	
%	 normUpdate = max([norm(dV)/norm(Vs(eqnNodes)),norm(currb-curra)/norm(curra),abs(currb(2)-currb(kx-2))/abs(currb(2))]);
%        normUpdate = max(norm(dV)/norm(Vs(eqnNodes)),norm(currb-curra)/norm(curra));

	normUpdate = max([norm(dV)/norm(Vs(eqnNodes)),norm(currb-curra)/norm(curra),abs(currt)/abs(currb(2))]);

	curra      = currb;

	display(['Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
	        '; ','normUpdate:',num2str(normUpdate),'; ','normdV:',num2str(norm(dV))]);
end  
% 
% Vs = Vs.*scl.Vt;
% ns = ns.*scl.ni;
% ps = ps.*scl.ni;

toc;   

display('End static solution cdbc');
