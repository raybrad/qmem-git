function [V,n,p,A] = simulSolutionAV(Vs,ns,ps,fren)

display(['Start dynamic solution at 1e',num2str(log10(fren)),' Hz']);

global sigma epsilon_in epsilon_sd;
global omega epsilon_mt;
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
V = zeros(Nnode,1);
V(dirNodes) = 1*acVolDirNodes;
n = zeros(Nnode,1);
p = zeros(Nnode,1);
A = zeros(Nlink,1);

eqnNodes = setdiff((1:Nnode)',dirNodes);
eqnSemiNodes = setdiff(semiNodes,[dirSemiNodes;dirNodes]);
eqnMetalNodes = setdiff(metalNodes,[bndNodes;dirSemiNodes]);
eqnLinks = setdiff((1:Nlink)',bndLinks);

%%%%%%%% start Newton's iteration %%%%%%%%
newtonTol = 1e-15;
updateTol = 1e-10;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 10;
maxLinIt = 200;
Eps = [epsilon_sd,epsilon_mt,epsilon_in];
Sgm = [0,sigma,0];
Eps_c = [epsilon_sd,sigma+1i*omega*epsilon_mt,epsilon_in];

NeqnNodes = length(eqnNodes);
NeqnSemiNodes = length(eqnSemiNodes);
Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
normRes_pre = 1e10;
itNr = 0;
tic;

curra = ones(kx,1);

while itNr < maxNewtonIt && normUpdate > updateTol
    
    JF_v = zeros(Nnode,Nnode); JF_n = zeros(Nnode,Nnode); 
    JF_p = zeros(Nnode,Nnode);
    JF_A = zeros(Nnode,Nlink);
    JKn_v = zeros(Nnode,Nnode); JKn_n = zeros(Nnode,Nnode); 
    JKn_p = zeros(Nnode,Nnode);
    JKn_A = zeros(Nnode,Nlink);
    JKp_v = zeros(Nnode,Nnode); JKp_n = zeros(Nnode,Nnode);
    JKp_p = zeros(Nnode,Nnode); JKp_A = zeros(Nnode,Nlink);
    JG_v = zeros(Nlink,Nnode); JG_n = zeros(Nlink,Nnode); 
    JG_p = zeros(Nlink,Nnode); 
    JG_A = zeros(Nlink,Nlink);
    rhs_F = zeros(Nnode,1); rhs_Kn = zeros(Nnode,1); 
    rhs_Kp = zeros(Nnode,1); 
    rhs_G = zeros(Nlink,1);
   
%    JF_p =[]; JKn_p =[]; JKp_v =[]; JKp_n =[]; JKp_p =[]; JKp_A =[]; JG_p =[];
%    rhs_Kp =[];

    currb=[];
    
    %%%%%%%% build Jacobian and rhs of F (Gauss's law) %%%%%%%%%%%%%
    for n1 = eqnNodes.'
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        sign_n1 = sign(ajnd_n1-n1);

            if any(ajvolM_n1 == 2) % metal/semi interfaces or triple points
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    ajvol_lk = linkVolumes{lk}(1,:);
                    ajvolS_lk = linkVolumes{lk}(2,:);
                    ajvolM_lk = volumeM(ajvol_lk);
                    E1 = -(V(n2)-V(n1))/linkL(lk)-sign_n1(i)*1i*omega*A(lk);
                    for j = 1:length(ajvol_lk)
                        switch ajvolM_lk(j)
                            case 1
                                alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                                beta = Vs(n2)-Vs(n1);
                                if isDirSemiNodes(n2), v2 = 0; else v2 = 1; end
                                sigma_n = sigma_c(mun,beta,'n',ns(n1),ns(n2));
                                sigma_p = sigma_c(mup,beta,'p',ps(n1),ps(n2));
                                
                                J_diff = Jc1_diff(alpha_n,beta,'n',n(n1),n(n2),0,v2)...
                                    +Jc1_diff(alpha_p,beta,'p',p(n1),p(n2),0,v2);
                                J_drift = (sigma_n+sigma_p)*E1;
                                Jd = 1i*omega*epsilon_sd*E1;
                                J = J_diff+J_drift+Jd;
                                dJv = (sigma_n+sigma_p+1i*omega*epsilon_sd)/linkL(lk);
                                dJA = -sign_n1(i)*1i*omega*(sigma_n+sigma_p+1i*omega*epsilon_sd);
                                dJn_nj = dJ1dc1j(alpha_n,beta,'n',v2);
                                JF_n(n1,n2) = JF_n(n1,n2)+ajvolS_lk(j)*dJn_nj;
                                dJp_pj = dJ1dc1j(alpha_p,beta,'p',v2);
                                JF_p(n1,n2) = JF_p(n1,n2)+ajvolS_lk(j)*dJp_pj;
                            case 2
                                J = (sigma+1i*omega*epsilon_mt)*E1;
                                dJv = (sigma+1i*omega*epsilon_mt)/linkL(lk);
                                dJA = -sign_n1(i)*1i*omega*(sigma+1i*omega*epsilon_mt);
                            case 3
                                J = 1i*omega*epsilon_in*E1;
                                dJv = 1i*omega*epsilon_in/linkL(lk);
                                dJA = sign_n1(i)*omega^2*epsilon_in;
                            otherwise
                                error('undefined material');
                        end
                        JF_v(n1,n1) = JF_v(n1,n1)+ajvolS_lk(j)*dJv;
                        JF_v(n1,n2) = JF_v(n1,n2)-ajvolS_lk(j)*dJv;
                        JF_A(n1,lk) = JF_A(n1,lk)+ajvolS_lk(j)*dJA;
                        rhs_F(n1) = rhs_F(n1)+ajvolS_lk(j)*J;
                    end
                end
            else
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    ajvol_lk = linkVolumes{lk}(1,:);
                    ajvolS_lk = linkVolumes{lk}(2,:);
                    ajvolM_lk = volumeM(ajvol_lk);
                    coefV = sum(Eps_c(ajvolM_lk).*ajvolS_lk)/linkL(lk);
                    JF_v(n1,n1) = JF_v(n1,n1)+coefV;
                    JF_v(n1,n2) = -coefV;
                    coefA = -1i*omega*coefV*linkL(lk);
                    JF_A(n1,lk) = sign_n1(i)*coefA;
                    E1 = -(V(n2)-V(n1))/linkL(lk)-sign_n1(i)*1i*omega*A(lk);
                    rhs_F(n1) = rhs_F(n1)+sum(Eps_c(ajvolM_lk).*ajvolS_lk)*E1;
                end
                semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
                JF_n(n1,n1) = semiV;
                JF_p(n1,n1) = -semiV;
                rhs_F(n1) = rhs_F(n1)+(n(n1)-p(n1))*semiV;          
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
            if isDirSemiNodes(n1), v1 = 0; else v1 = 1; end
            sign_n1 = sign(ajnd_n1-n1);
            for i = 1:length(ajlk_n1)
                n2 = ajnd_n1(i);
                lk = ajlk_n1(i);
                ajvol_lk = linkVolumes{lk}(1,:);
                ajvolS_lk = linkVolumes{lk}(2,:);
                ajvolM_lk = volumeM(ajvol_lk);
                semiS = sum(ajvolS_lk(ajvolM_lk == 1));
                alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                beta = Vs(n2)-Vs(n1);
                sigma_n = sigma_c(mun,beta,'n',ns(n1),ns(n2));
                sigma_p = sigma_c(mup,beta,'p',ps(n1),ps(n2));
                E1 = -(V(n2)-V(n1))/linkL(lk)-sign_n1(i)*1i*omega*A(lk);
                if isDirSemiNodes(n2), v2 = 0; else v2 = 1; end
                
                dJn_v = sigma_n/linkL(lk);
                JKn_v(n1,n1) = JKn_v(n1,n1)+semiS*dJn_v;
                JKn_v(n1,n2) = JKn_v(n1,n2)-semiS*dJn_v;
                dJp_v = sigma_p/linkL(lk);
                JKp_v(n1,n1) = JKp_v(n1,n1)+semiS*dJp_v;
                JKp_v(n1,n2) = JKp_v(n1,n2)-semiS*dJp_v;
                
                dJn_ni = dJ1dc1i(alpha_n,beta,'n',v1);
                dJn_nj = dJ1dc1j(alpha_n,beta,'n',v2);
                JKn_n(n1,n1) = JKn_n(n1,n1)+semiS*dJn_ni;
                JKn_n(n1,n2) = JKn_n(n1,n2)+semiS*dJn_nj;
                dJp_pi = dJ1dc1i(alpha_p,beta,'p',v1);
                dJp_pj = dJ1dc1j(alpha_p,beta,'p',v2);
                JKp_p(n1,n1) = JKp_p(n1,n1)+semiS*dJp_pi;
                JKp_p(n1,n2) = JKp_p(n1,n2)+semiS*dJp_pj;
                
                dJn_A = -1i*omega*sigma_n;
                dJp_A = -1i*omega*sigma_p;
                JKn_A(n1,lk) = sign_n1(i)*semiS*dJn_A;
                JKp_A(n1,lk) = sign_n1(i)*semiS*dJp_A;
                
                Jn = Jc1_diff(alpha_n,beta,'n',n(n1),n(n2),v1,v2)+sigma_n*E1;
                Jp = Jc1_diff(alpha_p,beta,'p',p(n1),p(n2),v1,v2)+sigma_p*E1;
                rhs_Kn(n1) = rhs_Kn(n1)+Jn*semiS;
                rhs_Kp(n1) = rhs_Kp(n1)+Jp*semiS;
            end
            semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
            JKn_n(n1,n1) = JKn_n(n1,n1)-1i*omega*semiV;
            JKp_p(n1,n1) = JKp_p(n1,n1)+1i*omega*semiV;
            rhs_Kn(n1) = rhs_Kn(n1)-1i*omega*semiV*n(n1);
            rhs_Kp(n1) = rhs_Kp(n1)+1i*omega*semiV*p(n1);
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
            coefV_n1 = scl.K*1i*omega*(sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1))*linkS(l1)/linkL(l1);
            JG_v(l1,n1) = JG_v(l1,n1)+coefV_n1;
            
            rhs_G(l1) = rhs_G(l1)-(-sign_n1.*coef_n1')*A(ajlk_n1);
            rhs_G(l1) = rhs_G(l1)+coefV_n1*V(n1);
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
            coefV_n2 = scl.K*1i*omega*(sum(Eps(ajvolM_n2).*ajvolV_n2)/nodeV(n2))*linkS(l1)/linkL(l1);
            JG_v(l1,n2) = JG_v(l1,n2)-coefV_n2;
            
            rhs_G(l1) = rhs_G(l1)-sign_n2.*coef_n2'*A(ajlk_n2);
            rhs_G(l1) = rhs_G(l1)-coefV_n2*V(n2);
        end
        
        JG_A(l1,:) = CC-GD;
        
        %%%%%% source current %%%%%%%%%%%
        E1 = -(V(n2)-V(n1))/linkL(l1)-1i*omega*A(l1);    
        for i = 1:length(ajvol_l1)
            switch ajvolM_l1(i)
                case 1
                    alpha_n = mun/linkL(l1); alpha_p = mup/linkL(l1);
                    beta = Vs(n2)-Vs(n1);
                    sigma_n = sigma_c(mun,beta,'n',ns(n1),ns(n2));
                    sigma_p = sigma_c(mup,beta,'p',ps(n1),ps(n2));
                    if isDirSemiNodes(n1), v1 = 0; else v1 = 1; end
                    if isDirSemiNodes(n2), v2 = 0; else v2 = 1; end
                    
                    J_diff = Jc1_diff(alpha_n,beta,'n',n(n1),n(n2),v1,v2)...
                        +Jc1_diff(alpha_p,beta,'p',p(n1),p(n2),v1,v2);
                    J_drift = (sigma_n+sigma_p)*E1;
                    J = J_diff+J_drift;
                    J = J+1i*omega*epsilon_sd*E1;
                    dJ_v = (sigma_n+sigma_p)/linkL(l1)+1i*omega*epsilon_sd/linkL(l1);
                    dJn_ni = dJ1dc1i(alpha_n,beta,'n',v1);
                    dJn_nj = dJ1dc1j(alpha_n,beta,'n',v2);
                    dJp_pi = dJ1dc1i(alpha_p,beta,'p',v1);
                    dJp_pj = dJ1dc1j(alpha_p,beta,'p',v2);
                    dJ_A = -1i*omega*(sigma_n+sigma_p)+omega^2*epsilon_sd;
                    JG_n(l1,n1) = JG_n(l1,n1)-scl.K*ajvolS_l1(i)*dJn_ni;
                    JG_n(l1,n2) = JG_n(l1,n2)-scl.K*ajvolS_l1(i)*dJn_nj;
                    JG_p(l1,n1) = JG_p(l1,n1)-scl.K*ajvolS_l1(i)*dJp_pi;
                    JG_p(l1,n2) = JG_p(l1,n2)-scl.K*ajvolS_l1(i)*dJp_pj;
                    
                case 2
                    dJ_v = (sigma+1i*omega*epsilon_mt)/linkL(l1);
                    dJ_A = -1i*omega*(sigma+1i*omega*epsilon_mt);
                    J = (sigma+1i*omega*epsilon_mt)*E1;
                case 3
                    dJ_v = 1i*omega*epsilon_in/linkL(l1);
                    dJ_A = omega^2*epsilon_in;
                    J = 1i*omega*epsilon_in*E1;
                otherwise
                    error('undefined material');
            end
            JG_v(l1,n1) = JG_v(l1,n1)-scl.K*ajvolS_l1(i)*dJ_v;
            JG_v(l1,n2) = JG_v(l1,n2)+scl.K*ajvolS_l1(i)*dJ_v;
            JG_A(l1,l1) = JG_A(l1,l1)-scl.K*ajvolS_l1(i)*dJ_A;
            rhs_G(l1) = rhs_G(l1)-scl.K*ajvolS_l1(i)*J;
        end
    end
    

    JF_v = JF_v(eqnNodes,eqnNodes);
    JF_n = JF_n(eqnNodes,eqnSemiNodes);
    JF_p = JF_p(eqnNodes,eqnSemiNodes);
    JF_A  = JF_A(eqnNodes,eqnLinks);
    JKn_v = JKn_v(eqnSemiNodes,eqnNodes);
    JKn_n = JKn_n(eqnSemiNodes,eqnSemiNodes);
    JKn_p = JKn_p(eqnSemiNodes,eqnSemiNodes);
    JKn_A  = JKn_A(eqnSemiNodes,eqnLinks);
    JKp_v = JKp_v(eqnSemiNodes,eqnNodes);
    JKp_n = JKp_n(eqnSemiNodes,eqnSemiNodes);
    JKp_p = JKp_p(eqnSemiNodes,eqnSemiNodes);
    JKp_A  = JKp_A(eqnSemiNodes,eqnLinks);
    JG_v = JG_v(eqnLinks,eqnNodes);
    JG_n = JG_n(eqnLinks,eqnSemiNodes);
    JG_p = JG_p(eqnLinks,eqnSemiNodes);
    JG_A  = JG_A(eqnLinks,eqnLinks);
    rhs_F = rhs_F(eqnNodes);
    rhs_Kn = rhs_Kn(eqnSemiNodes);
    rhs_Kp = rhs_Kp(eqnSemiNodes);
    rhs_G = rhs_G(eqnLinks);

    JF = [JF_v,JF_n,JF_p,JF_A];
    JKn = [JKn_v,JKn_n,JKn_p,JKn_A];
    JKp = [JKp_v,JKp_n,JKp_p,JKp_A];
    JG = [JG_v,JG_n,JG_p,JG_A];

    clear JF_v JF_n JF_p JF_A;
    clear JKn_v JKn_n JKn_p JKn_A;
    clear JKp_v JKp_n JKp_p JKp_A;
    clear JG_v JG_n JG_p JG_A;
      
    Jacob = [JF;JKn;JKp;JG];
    rhs = -[rhs_F;rhs_Kn;rhs_Kp;rhs_G];
    
    clear JF JKn JKp JG;
    clear rhs_F rhs_Kn rhs_Kp rhs_G;
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
%    dA = dX(NeqnNodes+NeqnSemiNodes+1:NeqnNodes+NeqnSemiNodes+Nl);
    normRes = norm(rhs);
%    if normRes > normRes_pre, error('Divergence detected'); end
    normRes_pre = normRes;

%    normUpdate = max([norm(dV)/norm(V(eqnNodes)),norm(dn)/norm(n(eqnSemiNodes)),...
%        norm(dp)/norm(p(eqnSemiNodes)),norm(dA)/norm(A(eqnLinks))]);
    
    

%    normUpdate = max([norm(dV)/norm(V(eqnNodes)),norm(dA)/norm(A(eqnLinks))]);

%    normUpdate = max([norm(dV)/norm(V(eqnNodes)), norm(dA)/norm(A(eqnLinks))]);    

    itNr = itNr+1;
    
    V(eqnNodes) = V(eqnNodes)+dV;
    n(eqnSemiNodes) = n(eqnSemiNodes)+dn;
    p(eqnSemiNodes) = p(eqnSemiNodes)+dp;
    A(eqnLinks) = A(eqnLinks)+dA;

    for i=1:kx
        currb=[currb;currentSav([i,5,5],Vs,ns,ps,V,n,p,A,1)];
    end

    currt=currentSav([1,round(ky/2),round(kz/2)],Vs,ns,ps,V,n,p,A,1)-currentSav([kx,round(ky/2),round(kz/2)],Vs,ns,ps,V,n,p,A,1)+currentSav([round(kx/2),1,round(kz/2)],Vs,ns,ps,V,n,p,A,2)-currentSav([round(kx/2),ky,round(kz/2)],Vs,ns,ps,V,n,p,A,2)+currentSav([round(kx/2),round(ky/2),1],Vs,ns,ps,V,n,p,A,3)-currentSav([round(kx/2),round(ky/2),kz],Vs,ns,ps,V,n,p,A,3);

    normUpdate = max([norm(dV)/norm(V(eqnNodes)),norm(currb-curra)/norm(curra),norm(dA)/norm(A(eqnLinks)),abs(currt)/abs(currb(2))]);


    curra      = currb;

    display(['Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
        '; ','normUpdate:',num2str(normUpdate),'; ','normdV:',num2str(norm(dV))]);

    clear dX dV dn dp dA;
    
end    
toc;   


%%%%%%%%%% postprocessing of variables %%%%%%%%%%%
% V_ = V.*scl.Vt;
% %%%% displaySlice(abs(V_),kz+1,'z');
% err_V = norm(V_-V_ref)/norm(V_ref);
% 
% if ~isempty(semiNodes)
%     n_ = n.*scl.ni;
%     err_n = norm(n_-n_ref)/norm(n_ref);
%     
%     p_ = p.*scl.ni;
%     err_p = norm(p_-p_ref)/norm(p_ref);
% end
% 
% 
% A_ = A.*scl.s_A;
% err_A = norm(A_-Aref_lk)/norm(Aref_lk);

display(['End dynamic solution at 1e',num2str(log10(fren)),' Hz']);

