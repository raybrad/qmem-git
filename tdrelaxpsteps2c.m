function [V,n,p,A,H,mJ,dtV,dtn,dtp,dtH] = tdrelaxpstep2(V_p,n_p,p_p,A_p,H_p,Js,mJ_p,dtV,dtn,dtp,dtH,dt)

display(['  Start TD simulation relaxation2c per step:']);

global sigma epsilon_in epsilon_sd;
global epsilon_mt;
global epsilon_qm;
global mun mup omega_p gamma_p;
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

global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global isqmnodes;

global currdlink;

global ehcrf;
global isqmvolm;

global light_speed XZsurfLinks XYsurfLinks YZsurfLinks;
global EsurfLinks BsurfLinks;
%%%%%%% initial guess  %%%%%%%%%%%%%
V = V_p + dtV *dt;
n = n_p + dtn *dt;
p = p_p + dtp *dt;
H = H_p + dtH *dt;
A = A_p + H_p *dt + dtH * dt^2/2;
mJ=mJ_p;

dtV1 = dtV;
dtn1 = dtn;
dtp1 = dtp;
dtH1 = dtH;

eqnNodes = setdiff((1:Nnode)',[dirNodes]);
eqnSemiNodes = setdiff(semiNodes,[dirSemiNodes;dirNodes;sQMnodes]);

eqnLinks=setdiff((1:Nlink)',EsurfLinks);     %for E boundary condition

n(sQMnodes)=n_p(sQMnodes);
p(sQMnodes)=p_p(sQMnodes);

isBsurfLinks=false(Nlink,1);
isBsurfLinks(BsurfLinks)=true;		     %for B boundary condition
%%%%%%%% start Newton's iteration %%%%%%%%
newtonTol = 1e-15;
updateTol = 1e-10;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 10;
maxLinIt = 200;
Eps = [epsilon_sd,epsilon_mt,epsilon_in,epsilon_qm];
Sgm = [0,sigma,0];

NeqnNodes = length(eqnNodes);
NeqnSemiNodes = length(eqnSemiNodes);
Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
normRes_pre = 1e10;
itNr = 0;
tStart=tic;
while itNr < maxNewtonIt && normUpdate > updateTol

     rhs_F = zeros(Nnode,1); rhs_Kn = zeros(Nnode,1);
     rhs_Kp =zeros(Nnode,1);
     rhs_G = zeros(Nlink,1);

     ntripletsJFv=48*Nnode;
     rowJFv=zeros(ntripletsJFv,1);
     colJFv=zeros(ntripletsJFv,1);
     valJFv=zeros(ntripletsJFv,1);
     ntripletsJFv=0;

     ntripletsJFH=24*Nnode;
     rowJFH=zeros(ntripletsJFH,1);
     colJFH=zeros(ntripletsJFH,1);
     valJFH=zeros(ntripletsJFH,1);
     ntripletsJFH=0;

     ntripletsJFn=24*Nnode;
     rowJFn=zeros(ntripletsJFn,1);
     colJFn=zeros(ntripletsJFn,1);
     valJFn=zeros(ntripletsJFn,1);
     ntripletsJFn=0;

     ntripletsJFp=24*Nnode;
     rowJFp=zeros(ntripletsJFp,1);
     colJFp=zeros(ntripletsJFp,1);
     valJFp=zeros(ntripletsJFp,1);
     ntripletsJFp=0;

    %update mJ,for plasmonic metal
    %mJ_(n+1)=(2-gamma_p*dt)/(2+gamma_p*dt)*mJ_n+ omega_p^2*dt*(E_n+1 + E_n)/(2+gamma_p*dt)
    for l1 = eqnLinks'
    n1 = links(l1,1);
    n2 = links(l1,2);
    ajvol_l1 = linkVolumes{l1}(1,:);%volumeid
    ajvolS_l1 = linkVolumes{l1}(2,:);%associate surf area,not dual area, but only its own 1/4
    ajvolM_l1 = volumeM(ajvol_l1);
    for i = 1:length(ajvol_l1)
     if  (ajvolM_l1(i)==2)
     mJ(l1) =(2-gamma_p*dt)/(2+gamma_p*dt)*mJ_p(l1)+dt*omega_p*omega_p*((V(n1)+V_p(n1)-V(n2)-V_p(n2))/linkL(l1)-H(l1)-H_p(l1))/(2+gamma_p*dt);
     end
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     for k=1:NeqnNodes
        n1 = eqnNodes(k);
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        for qmi = 1:length(ajvol_n1)
          if isqmvolm(ajvol_n1(qmi))
           ajvolM_n1(qmi)=4;	%
          end 
        end
        sign_n1 = sign(ajnd_n1-n1);

        if isDirSemiNodes(n1), dV1 = deltaV(doping(n1)); else dV1 = 0; end

        if all(ajvolM_n1 == 4)	%all nodes in QM region (include the links along the boundary)
            for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    n3 = links(lk,3);
                    dtE1 = -(dtV(n2)-dtV(n1))/linkL(lk)-sign_n1(i)*dtH(lk);
                    rhs_F(n1)   = rhs_F(n1) + sign_n1(i)*linkS(lk)*currdlink(lk,n3+3) + linkS(lk)* epsilon_qm *dtE1;
		    ntripletsJFv=ntripletsJFv+1;
		    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n1;
		    valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+epsilon_qm *linkS(lk)/(0.5*dt*linkL(lk));
		    ntripletsJFv=ntripletsJFv+1;
		    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n2;
		    valJFv(ntripletsJFv)=valJFv(ntripletsJFv)-epsilon_qm *linkS(lk)/(0.5*dt*linkL(lk));
		    ntripletsJFH=ntripletsJFH+1;
		    rowJFH(ntripletsJFH)=n1;
		    colJFH(ntripletsJFH)=lk;
		    valJFH(ntripletsJFH)=valJFH(ntripletsJFH)-epsilon_qm *linkS(lk)*sign_n1(i)/(0.5*dt);
            end
        elseif any(ajvolM_n1 == 2) % metal/semi interfaces or triple points
            for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    n3 = links(lk,3);
                    dtE1 = -(dtV(n2)-dtV(n1))/linkL(lk)-sign_n1(i)*dtH(lk);
                    if isSQMlinks(lk) == 0	%along(impossible here) or connected to the outside QM boundary,or pure EM links 
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

		    ntripletsJFn=ntripletsJFn+1;
		    ntripletsJFp=ntripletsJFp+1;
		    rowJFn(ntripletsJFn)=n1;
		    colJFn(ntripletsJFn)=n2;
		    rowJFp(ntripletsJFp)=n1;
		    colJFp(ntripletsJFp)=n2;
                      for j = 1:length(ajvol_lk)
                         switch ajvolM_lk(j)
                            case 1
                                alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                                beta = V(n2)-V(n1) + sign_n1(i)*H(lk)*linkL(lk);
                                J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                                    +Jc(alpha_p,beta,'p',p(n1),p(n2));
                                if isqmvolm(ajvol_lk(j))
                                  epsilon_cu = epsilon_qm;
                                else 
                                  epsilon_cu = epsilon_sd;
                                end
                                Jd = epsilon_cu * dtE1;
                                J = J_diff +Jd;
                                dJv = dJndV(alpha_n,beta,n(n1),n(n2))...
                                    +dJpdV(alpha_p,beta,p(n1),p(n2)) + epsilon_cu/(0.5*dt*linkL(lk));
                                dJH = -sign_n1(i)*dJv*linkL(lk);
                                dJn_nj = dJdcj(alpha_n,beta,'n');
				valJFn(ntripletsJFn)=valJFn(ntripletsJFn)+ttt*ajvolS_lk(j)*dJn_nj;
                                dJp_pj = dJdcj(alpha_p,beta,'p');
				valJFp(ntripletsJFp)=valJFp(ntripletsJFp)+ttt*ajvolS_lk(j)*dJp_pj;
                            case 2
                  		J=sign_n1(i)*mJ(lk)+epsilon_mt*dtE1;
                                dJv = sign_n1(i)*dt*omega_p*omega_p/(dt*gamma_p+2)/linkL(lk)+ epsilon_mt/(0.5*dt*linkL(lk));
                                dJH = -sign_n1(i)*(dt*omega_p*omega_p/(dt*gamma_p+2)+epsilon_mt/(0.5*dt));
                            case 3
                                J   = epsilon_in*dtE1;
                                dJv = epsilon_in/(0.5*dt*linkL(lk));
                                dJH = -sign_n1(i)*epsilon_in/(0.5*dt);
                            otherwise
                                error('undefined material');
                        end
			ntripletsJFv=ntripletsJFv+1;
			rowJFv(ntripletsJFv)=n1;
			colJFv(ntripletsJFv)=n1;
			valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+ajvolS_lk(j)*dJv;
			ntripletsJFv=ntripletsJFv+1;
			rowJFv(ntripletsJFv)=n1;
			colJFv(ntripletsJFv)=n2;
			valJFv(ntripletsJFv)=valJFv(ntripletsJFv)-ajvolS_lk(j)*dJv;

			ntripletsJFH=ntripletsJFH+1;
			rowJFH(ntripletsJFH)=n1;
			colJFH(ntripletsJFH)=lk;
			valJFH(ntripletsJFH)=valJFH(ntripletsJFH)+ajvolS_lk(j)*dJH;

                        rhs_F(n1) = rhs_F(n1)+ajvolS_lk(j)*J;
                      end
                    elseif isSQMlinks(lk) <= 6	%inside/outside interface
                       rhs_F(n1)   = rhs_F(n1) + sign_n1(i)*linkS(lk)*currdlink(lk,n3+3) + linkS(lk)* epsilon_qm *dtE1;
		       ntripletsJFv=ntripletsJFv+1;
		       rowJFv(ntripletsJFv)=n1;
		       colJFv(ntripletsJFv)=n1;
		       valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+ epsilon_qm *linkS(lk)/(0.5*dt*linkL(lk)); 
		       ntripletsJFv=ntripletsJFv+1;
		       rowJFv(ntripletsJFv)=n1;
		       colJFv(ntripletsJFv)=n2;
		       valJFv(ntripletsJFv)=valJFv(ntripletsJFv)- epsilon_qm *linkS(lk)/(0.5*dt*linkL(lk)); 

		       ntripletsJFH=ntripletsJFH+1;
		       rowJFH(ntripletsJFH)=n1;
		       colJFH(ntripletsJFH)=lk;
		       valJFH(ntripletsJFH)=valJFH(ntripletsJFH)- epsilon_qm *linkS(lk)*sign_n1(i)/(0.5*dt);
                    elseif isSQMlinks(lk) > 6 %inside QM, very rare in ajvolM_n1 == 2 not ajvolM_n1==4 ,it's 0 anyway
                       rhs_F(n1)   = rhs_F(n1) + 0;
                    end
            end
        else   %other ajvolM_n1 except QMnodes and metal 
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    if isSQMlinks(lk) == 0  %along(impossible here)or connected to the outside QM boundary,or pure EM links 
                      if isDirSemiNodes(n2)
                        dV2 = deltaV(doping(n2));
                      else
                        dV2 = 0;
                      end
                      ajvol_lk = linkVolumes{lk}(1,:);
                      ajvolS_lk = linkVolumes{lk}(2,:);
                      ajvolM_lk = volumeM(ajvol_lk);
                      for qmi = 1:length(ajvol_lk)
                         if isqmvolm(ajvol_lk(qmi))
                          ajvolM_lk(qmi)=4;
                         end
                      end
                      coefV = sum(Eps(ajvolM_lk).*ajvolS_lk)/linkL(lk);
		      ntripletsJFv=ntripletsJFv+1;
		      rowJFv(ntripletsJFv)=n1;
		      colJFv(ntripletsJFv)=n1;
		      valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+coefV;
		      ntripletsJFv=ntripletsJFv+1;
		      rowJFv(ntripletsJFv)=n1;
		      colJFv(ntripletsJFv)=n2;
		      valJFv(ntripletsJFv)=-coefV;
                      dV2_vec = [0,1,1/2,0];
                      coefH = -coefV*linkL(lk);
		      ntripletsJFH=ntripletsJFH+1;
		      rowJFH(ntripletsJFH)=n1;
		      colJFH(ntripletsJFH)=lk;
		      valJFH(ntripletsJFH)=sign_n1(i)*coefH;
                      rhs_F(n1) = rhs_F(n1)+sum(Eps(ajvolM_lk).*(-(V(n2)-V(n1)+dV2_vec(ajvolM_lk).*dV2)/linkL(lk)-sign_n1(i)* H(lk)).*ajvolS_lk);
                    else	%interface
                     rhs_F(n1)   = rhs_F(n1) + epsilon_qm *(-(V(n2)-V(n1))/linkL(lk)-sign_n1(i)* H(lk))*linkS(lk);
		     ntripletsJFv=ntripletsJFv+1;
		     rowJFv(ntripletsJFv)=n1;
		     colJFv(ntripletsJFv)=n1;
		     valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+epsilon_qm * linkS(lk)/linkL(lk);
		     ntripletsJFv=ntripletsJFv+1;
		     rowJFv(ntripletsJFv)=n1;
		     colJFv(ntripletsJFv)=n2;
		     valJFv(ntripletsJFv)=valJFv(ntripletsJFv)-epsilon_qm * linkS(lk)/linkL(lk);
		     ntripletsJFH=ntripletsJFH+1;
		     rowJFH(ntripletsJFH)=n1;
		     colJFH(ntripletsJFH)=lk;
		     valJFH(ntripletsJFH)=  -sign_n1(i) * epsilon_qm *linkS(lk); 
                    end
                end
                semiV = sum(ajvolV_n1(ajvolM_n1 == 1)) + sum(ajvolV_n1(ajvolM_n1 == 4));
		ntripletsJFn=ntripletsJFn+1;
		rowJFn(ntripletsJFn)=n1;
		colJFn(ntripletsJFn)=n1;
		valJFn(ntripletsJFn)=semiV;
		ntripletsJFp=ntripletsJFp+1;
		rowJFp(ntripletsJFp)=n1;
		colJFp(ntripletsJFp)=n1;
		valJFp(ntripletsJFp)=-semiV;
                rhs_F(n1) = rhs_F(n1)+(n(n1)-p(n1) - doping(n1))*semiV;
        end
%        if any(ajvolM_n1 == 2) % scaling some equations with condutivity for matrix balance
%            JF_v(n1,:) = JF_v(n1,:)/sigma;
%            JF_n(n1,:) = JF_n(n1,:)/sigma;
%            JF_p(n1,:) = JF_p(n1,:)/sigma;
%            JF_H(n1,:) = JF_H(n1,:)/sigma;
%            rhs_F(n1) = rhs_F(n1)/sigma;
%        end
     end
JF_v=sparse(rowJFv(1:ntripletsJFv),colJFv(1:ntripletsJFv),valJFv(1:ntripletsJFv),Nnode,Nnode);
JF_H=sparse(rowJFH(1:ntripletsJFH),colJFH(1:ntripletsJFH),valJFH(1:ntripletsJFH),Nnode,Nlink);
JF_n=sparse(rowJFn(1:ntripletsJFn),colJFn(1:ntripletsJFn),valJFn(1:ntripletsJFn),Nnode,Nnode);
JF_p=sparse(rowJFp(1:ntripletsJFp),colJFp(1:ntripletsJFp),valJFp(1:ntripletsJFp),Nnode,Nnode);
clear rowJFv colJFv valJFv rowJFH colJFH valJFH rowJFn colJFn valJFn rowJFp colJFp valJFp;
JF_v = JF_v(eqnNodes,eqnNodes);
JF_n = JF_n(eqnNodes,eqnSemiNodes);
JF_p = JF_p(eqnNodes,eqnSemiNodes);
JF_H = JF_H(eqnNodes,eqnLinks);
JF = [JF_v,JF_n,JF_p,JF_H];
%save('dump/JF.mat','JF_v','JF_n','JF_p','JF_H');
clear JF_v JF_n JF_p JF_H;
display(['time for matrix collection,F matrix:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%% build Jacobian and rhs of K (current continuity of semi) %%%%%%%%%%%%%
     if ~isempty(eqnSemiNodes)
     ntripletsJKnv=12*NeqnSemiNodes;
     rowJKnv=zeros(ntripletsJKnv,1);
     colJKnv=zeros(ntripletsJKnv,1);
     valJKnv=zeros(ntripletsJKnv,1);
     ntripletsJKnv=0;

     ntripletsJKnH=6*NeqnSemiNodes;
     rowJKnH=zeros(ntripletsJKnH,1);
     colJKnH=zeros(ntripletsJKnH,1);
     valJKnH=zeros(ntripletsJKnH,1);
     ntripletsJKnH=0;

     ntripletsJKnn=12*NeqnSemiNodes;
     rowJKnn=zeros(ntripletsJKnn,1);
     colJKnn=zeros(ntripletsJKnn,1);
     valJKnn=zeros(ntripletsJKnn,1);
     ntripletsJKnn=0;

     ntripletsJKpv=12*NeqnSemiNodes;
     rowJKpv=zeros(ntripletsJKpv,1);
     colJKpv=zeros(ntripletsJKpv,1);
     valJKpv=zeros(ntripletsJKpv,1);
     ntripletsJKpv=0;

     ntripletsJKpH=6*NeqnSemiNodes;
     rowJKpH=zeros(ntripletsJKpH,1);
     colJKpH=zeros(ntripletsJKpH,1);
     valJKpH=zeros(ntripletsJKpH,1);
     ntripletsJKpH=0;

     ntripletsJKpp=12*NeqnSemiNodes;
     rowJKpp=zeros(ntripletsJKpp,1);
     colJKpp=zeros(ntripletsJKpp,1);
     valJKpp=zeros(ntripletsJKpp,1);
     ntripletsJKpp=0;
        for k=1:NeqnSemiNodes
	    n1 = eqnSemiNodes(k);
            ajlk_n1 = nodeLinks{n1}(1,:);
            ajnd_n1 = nodeLinks{n1}(2,:);
            ajvol_n1 = nodeVolumes{n1}(1,:);
            ajvolV_n1 = nodeVolumes{n1}(2,:);
            ajvolM_n1 = volumeM(ajvol_n1);
            sign_n1 = sign(ajnd_n1-n1);
            for i = 1:length(ajlk_n1)
                n2 = ajnd_n1(i);
                lk = ajlk_n1(i);
                if isSQMlinks(lk) == 0  %along or connected to the outside QM boundary,or pure EM links(impossible here)
                  ajvol_lk = linkVolumes{lk}(1,:);
                  ajvolS_lk = linkVolumes{lk}(2,:);
                  ajvolM_lk = volumeM(ajvol_lk);
                  semiS = sum(ajvolS_lk(ajvolM_lk == 1));
                  alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                  E1 = -(V(n2)-V(n1))/linkL(lk)-sign_n1(i)*H(lk);
                  beta = - linkL(lk) * E1;

                  dJn_v = dJndV(alpha_n,beta,n(n1),n(n2));
	    	  ntripletsJKnv=ntripletsJKnv+1;
		  rowJKnv(ntripletsJKnv)=n1;
		  colJKnv(ntripletsJKnv)=n1;
		  valJKnv(ntripletsJKnv)=valJKnv(ntripletsJKnv)+semiS*dJn_v;
	    	  ntripletsJKnv=ntripletsJKnv+1;
		  rowJKnv(ntripletsJKnv)=n1;
		  colJKnv(ntripletsJKnv)=n2;
		  valJKnv(ntripletsJKnv)=valJKnv(ntripletsJKnv)-semiS*dJn_v;

                  dJp_v = dJpdV(alpha_p,beta,p(n1),p(n2));
	    	  ntripletsJKpv=ntripletsJKpv+1;
		  rowJKpv(ntripletsJKpv)=n1;
		  colJKpv(ntripletsJKpv)=n1;
		  valJKpv(ntripletsJKpv)=valJKpv(ntripletsJKpv)+semiS*dJp_v;
	    	  ntripletsJKpv=ntripletsJKpv+1;
		  rowJKpv(ntripletsJKpv)=n1;
		  colJKpv(ntripletsJKpv)=n2;
		  valJKpv(ntripletsJKpv)=valJKpv(ntripletsJKpv)-semiS*dJp_v;

                  dJn_ni = dJdci(alpha_n,beta,'n');
                  dJn_nj = dJdcj(alpha_n,beta,'n');
	    	  ntripletsJKnn=ntripletsJKnn+1;
		  rowJKnn(ntripletsJKnn)=n1;
		  colJKnn(ntripletsJKnn)=n1;
		  valJKnn(ntripletsJKnn)=valJKnn(ntripletsJKnn)+semiS*dJn_ni;
	    	  ntripletsJKnn=ntripletsJKnn+1;
		  rowJKnn(ntripletsJKnn)=n1;
		  colJKnn(ntripletsJKnn)=n2;
		  valJKnn(ntripletsJKnn)=valJKnn(ntripletsJKnn)+semiS*dJn_nj;

                  dJp_pi = dJdci(alpha_p,beta,'p');
                  dJp_pj = dJdcj(alpha_p,beta,'p');
	    	  ntripletsJKpp=ntripletsJKpp+1;
		  rowJKpp(ntripletsJKpp)=n1;
		  colJKpp(ntripletsJKpp)=n1;
		  valJKpp(ntripletsJKpp)=valJKpp(ntripletsJKpp)+semiS*dJp_pi;
	    	  ntripletsJKpp=ntripletsJKpp+1;
		  rowJKpp(ntripletsJKpp)=n1;
		  colJKpp(ntripletsJKpp)=n2;
		  valJKpp(ntripletsJKpp)=valJKpp(ntripletsJKpp)+semiS*dJp_pj;

	    	  ntripletsJKnH=ntripletsJKnH+1;
		  rowJKnH(ntripletsJKnH)=n1;
		  colJKnH(ntripletsJKnH)=lk;
		  valJKnH(ntripletsJKnH)=valJKnH(ntripletsJKnH)-sign_n1(i)*semiS*dJn_v*linkL(lk);

	    	  ntripletsJKpH=ntripletsJKpH+1;
		  rowJKpH(ntripletsJKpH)=n1;
		  colJKpH(ntripletsJKpH)=lk;
		  valJKpH(ntripletsJKpH)=valJKpH(ntripletsJKpH)-sign_n1(i)*semiS*dJp_v*linkL(lk);

                  Jn = Jc(alpha_n,beta,'n',n(n1),n(n2));
                  Jp = Jc(alpha_p,beta,'p',p(n1),p(n2));
                  rhs_Kn(n1) = rhs_Kn(n1)+Jn*semiS;
                  rhs_Kp(n1) = rhs_Kp(n1)+Jp*semiS;
                elseif isSQMlinks(lk) == 1
                    rhs_Kn(n1)   = rhs_Kn(n1) + ehcrf     * linkS(lk)*currdlink(lk,4);
                    rhs_Kp(n1)   = rhs_Kp(n1) + (1-ehcrf) * linkS(lk)*currdlink(lk,4);
                elseif isSQMlinks(lk) == 2
                    rhs_Kn(n1)   = rhs_Kn(n1) - ehcrf     * linkS(lk)*currdlink(lk,4);
                    rhs_Kp(n1)   = rhs_Kp(n1) - (1-ehcrf) * linkS(lk)*currdlink(lk,4);
                elseif isSQMlinks(lk) == 3
                    rhs_Kn(n1)   = rhs_Kn(n1) + ehcrf     * linkS(lk)*currdlink(lk,5);
                    rhs_Kp(n1)   = rhs_Kp(n1) + (1-ehcrf) * linkS(lk)*currdlink(lk,5);
                elseif isSQMlinks(lk) == 4
                    rhs_Kn(n1)   = rhs_Kn(n1) - ehcrf     * linkS(lk)*currdlink(lk,5);
                    rhs_Kp(n1)   = rhs_Kp(n1) - (1-ehcrf) * linkS(lk)*currdlink(lk,5);
                elseif isSQMlinks(lk) == 5
                    rhs_Kn(n1)   = rhs_Kn(n1) + ehcrf     * linkS(lk)*currdlink(lk,6);
                    rhs_Kp(n1)   = rhs_Kp(n1) + (1-ehcrf) * linkS(lk)*currdlink(lk,6);
                elseif isSQMlinks(lk) == 6
                    rhs_Kn(n1)   = rhs_Kn(n1) - ehcrf     * linkS(lk)*currdlink(lk,6);
                    rhs_Kp(n1)   = rhs_Kp(n1) - (1-ehcrf) * linkS(lk)*currdlink(lk,6);
                elseif isSQMlinks(lk) > 6
                    rhs_Kn(n1)   = rhs_Kn(n1) + 0;
                    rhs_Kp(n1)   = rhs_Kp(n1) + 0;
                end
            end
            semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
	    ntripletsJKnn=ntripletsJKnn+1;
	    rowJKnn(ntripletsJKnn)=n1;
	    colJKnn(ntripletsJKnn)=n1;
	    valJKnn(ntripletsJKnn)=valJKnn(ntripletsJKnn)-semiV/(0.5*dt);
	    ntripletsJKpp=ntripletsJKpp+1;
	    rowJKpp(ntripletsJKpp)=n1;
	    colJKpp(ntripletsJKpp)=n1;
	    valJKpp(ntripletsJKpp)=valJKpp(ntripletsJKpp)+semiV/(0.5*dt);
            rhs_Kn(n1) = rhs_Kn(n1)-semiV*dtn(n1);
            rhs_Kp(n1) = rhs_Kp(n1)+semiV*dtp(n1);
        end
     end

if ~isempty(eqnSemiNodes)
JKn_v=sparse(rowJKnv(1:ntripletsJKnv),colJKnv(1:ntripletsJKnv),valJKnv(1:ntripletsJKnv),Nnode,Nnode);
JKn_H=sparse(rowJKnH(1:ntripletsJKnH),colJKnH(1:ntripletsJKnH),valJKnH(1:ntripletsJKnH),Nnode,Nlink);
JKn_n=sparse(rowJKnn(1:ntripletsJKnn),colJKnn(1:ntripletsJKnn),valJKnn(1:ntripletsJKnn),Nnode,Nnode);
JKp_v=sparse(rowJKpv(1:ntripletsJKpv),colJKpv(1:ntripletsJKpv),valJKpv(1:ntripletsJKpv),Nnode,Nnode);
JKp_H=sparse(rowJKpH(1:ntripletsJKpH),colJKpH(1:ntripletsJKpH),valJKpH(1:ntripletsJKpH),Nnode,Nlink);
JKp_p=sparse(rowJKpp(1:ntripletsJKpp),colJKpp(1:ntripletsJKpp),valJKpp(1:ntripletsJKpp),Nnode,Nnode);
clear rowJKnv colJKnv valJKnv rowJKnH colJKnH valJKnH rowJKnn colJKnn valJKnn;
clear rowJKpv colJKpv valJKpv rowJKpH colJKpH valJKpH rowJKpp colJKpp valJKpp;
JKn_v = JKn_v(eqnSemiNodes,eqnNodes);
JKn_n = JKn_n(eqnSemiNodes,eqnSemiNodes);
JKn_p = [];
JKn_H  = JKn_H(eqnSemiNodes,eqnLinks);
JKp_v = JKp_v(eqnSemiNodes,eqnNodes);
JKp_p = JKp_p(eqnSemiNodes,eqnSemiNodes);
JKp_n = [];
JKp_H  = JKp_H(eqnSemiNodes,eqnLinks);
else
JKn_v=[];JKn_H=[];JKn_p=[];JKn_n=[];
JKp_v=[];JKp_H=[];JKp_n=[];JKp_p=[];
end
JKn = [JKn_v,JKn_n,JKn_p,JKn_H];
JKp = [JKp_v,JKp_n,JKp_p,JKp_H];
clear JKn_v JKn_n JKn_p JKn_H JKp_v JKp_p JKp_n JKp_H;
display(['time for matrix collection,K matrix:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
   ntripletsJGv=24*Nlink;     
   rowJGv=zeros(ntripletsJGv,1);  
   colJGv=zeros(ntripletsJGv,1);
   valJGv=zeros(ntripletsJGv,1);  
   ntripletsJGv=0;                
   ntripletsJGH=(12+4*4+30)*Nlink;      
   rowJGH=zeros(ntripletsJGH,1);  
   colJGH=zeros(ntripletsJGH,1);  
   valJGH=zeros(ntripletsJGH,1);  
   ntripletsJGH=0;            
   ntripletsJGn=24*Nlink;     
   rowJGn=zeros(ntripletsJGn,1);  
   colJGn=zeros(ntripletsJGn,1);  
   valJGn=zeros(ntripletsJGn,1);  
   ntripletsJGn=0;                
   ntripletsJGp=24*Nlink;      
   rowJGp=zeros(ntripletsJGp,1);  
   colJGp=zeros(ntripletsJGp,1);  
   valJGp=zeros(ntripletsJGp,1);  
   ntripletsJGp=0;               

    for k=1:Nl 
        l1 = eqnLinks(k);
        n1 = links(l1,1);
        n2 = links(l1,2);
        ajvol_l1 = linkVolumes{l1}(1,:);
        ajvolS_l1 = linkVolumes{l1}(2,:);
        ajvolM_l1 = volumeM(ajvol_l1);
        for qmi = 1:length(ajvol_l1)
          if isqmvolm(ajvol_l1(qmi))
           ajvolM_l1(qmi)=4;
          end
        end
        n3 = links(l1,3);
        CC = zeros(1,Nlink);
        GD = zeros(1,Nlink);
        bndJG_H = zeros(1,Nlink);
	ajlk_n1=[];
	ajlk_n2=[];
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
	    for m=1:length(ajlk)
	    ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=ajlk(m);
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+CC(ajlk(m))*dt/2;
	    end
        end

        rhs_G(l1) = rhs_G(l1)+CC(l1)*A(l1);
	ntripletsJGH=ntripletsJGH+1;
	rowJGH(ntripletsJGH)=l1;
	colJGH(ntripletsJGH)=l1;
	valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+CC(l1)*dt/2;

	    if (size(ajsf_l1,2)==3) %xy plane : Ay(lk)~Bx~par Ay(lk)/par z =1/c par Ay(lk)/par t no matter top or bottom
	        if(isBsurfLinks(l1))
		rhs_G(l1)=rhs_G(l1)+0.0;
		bndJG_H(l1)=0.0;
	        ntripletsJGH=ntripletsJGH+1;
		rowJGH(ntripletsJGH)=l1;
		colJGH(ntripletsJGH)=l1;
		valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
       	    else
       		rhs_G(l1)=rhs_G(l1)+dL*(1/light_speed*H(l1));
       		bndJG_H(l1)=dL/light_speed;
	        ntripletsJGH=ntripletsJGH+1;
		rowJGH(ntripletsJGH)=l1;
		colJGH(ntripletsJGH)=l1;
		valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
	    end
	    elseif (size(ajsf_l1,2)==2) 
	        if(isBsurfLinks(l1))
		rhs_G(l1)=rhs_G(l1)+dL*(1/light_speed*H(l1));
		bndJG_H(l1)=dL/light_speed;
	        ntripletsJGH=ntripletsJGH+1;
		rowJGH(ntripletsJGH)=l1;
		colJGH(ntripletsJGH)=l1;
		valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
		else
		rhs_G(l1)=rhs_G(l1)+2*dL*(1/light_speed*H(l1));
		bndJG_H(l1)=2*dL/light_speed;
	        ntripletsJGH=ntripletsJGH+1;
		rowJGH(ntripletsJGH)=l1;
		colJGH(ntripletsJGH)=l1;
		valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
		end
	    end
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
            for qmi = 1:length(ajvol_n1)
              if isqmvolm(ajvol_n1(qmi))
                ajvolM_n1(qmi)=4;
              end 
            end
            coefV_n1 = scl.K*(sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1))*linkS(l1)/linkL(l1);

            ntripletsJGv=ntripletsJGv+1;
	    rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n1;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+coefV_n1/(0.5*dt);

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
            for qmi = 1:length(ajvol_n2)
              if isqmvolm(ajvol_n2(qmi))
                ajvolM_n2(qmi)=4;
              end
            end
            coefV_n2 = scl.K*(sum(Eps(ajvolM_n2).*ajvolV_n2)/nodeV(n2))*linkS(l1)/linkL(l1);

            ntripletsJGv=ntripletsJGv+1;
	    rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n2;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)-coefV_n2/(0.5*dt);

            rhs_G(l1) = rhs_G(l1)-sign_n2.*coef_n2'*A(ajlk_n2);
            rhs_G(l1) = rhs_G(l1)-coefV_n2*dtV(n2);
        end

	ajlk_n1n2=unique([ajlk_n1;ajlk_n2]);
        for m=1:length(ajlk_n1n2)
	ntripletsJGH=ntripletsJGH+1;
	rowJGH(ntripletsJGH)=l1;
	colJGH(ntripletsJGH)=ajlk_n1n2(m);
	valJGH(ntripletsJGH)=valJGH(ntripletsJGH)-GD(ajlk_n1n2(m))*dt/2;
	end

        %%%%%% source current %%%%%%%%%%%
        E1 = -(V(n2)-V(n1))/linkL(l1)-H(l1); %(div V+ PI)
        dtE1 = -(dtV(n2)-dtV(n1))/linkL(l1) - dtH(l1);
        if isSQMlinks(l1) == 0		 %along or connected to the outside QM boundary,or pure EM links
          E1 = -(V(n2)-V(n1))/linkL(l1)-H(l1);
          if isDirSemiNodes(n1), dV1 = deltaV(doping(n1)); else dV1 = 0; end
          if isDirSemiNodes(n2), dV2 = deltaV(doping(n2)); else dV2 = 0; end
          for i = 1:length(ajvol_l1)	%volumetype of the link belong to
            switch ajvolM_l1(i)
                case 1
                    alpha_n = mun/linkL(l1); alpha_p = mup/linkL(l1);
                    beta = V(n2)-V(n1)+H(l1)*linkL(l1);
                    J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                        +Jc(alpha_p,beta,'p',p(n1),p(n2));
                    J = J_diff + epsilon_sd*dtE1;
                    dJ_v = dJndV(alpha_n,beta,n(n1),n(n2))...
                          +dJpdV(alpha_p,beta,p(n1),p(n2))+epsilon_sd/(0.5*dt*linkL(l1));
                    dJn_ni = dJdci(alpha_n,beta,'n');
                    dJn_nj = dJdcj(alpha_n,beta,'n');
                    dJp_pi = dJdci(alpha_p,beta,'p');
                    dJp_pj = dJdcj(alpha_p,beta,'p');
                    dJ_H = - dJ_v*linkL(l1);

        	    ntripletsJGn=ntripletsJGn+1;
              	    rowJGn(ntripletsJGn)=l1;
		    colJGn(ntripletsJGn)=n1;
		    valJGn(ntripletsJGn)=valJGn(ntripletsJGn)-scl.K*ajvolS_l1(i)*dJn_ni;
		    ntripletsJGn=ntripletsJGn+1;
		    rowJGn(ntripletsJGn)=l1;
		    colJGn(ntripletsJGn)=n2;
		    valJGn(ntripletsJGn)=valJGn(ntripletsJGn)-scl.K*ajvolS_l1(i)*dJn_nj;

        	    ntripletsJGp=ntripletsJGp+1;
              	    rowJGp(ntripletsJGp)=l1;
		    colJGp(ntripletsJGp)=n1;
		    valJGp(ntripletsJGp)=valJGp(ntripletsJGp)-scl.K*ajvolS_l1(i)*dJp_pi;
		    ntripletsJGp=ntripletsJGp+1;
		    rowJGp(ntripletsJGp)=l1;
		    colJGp(ntripletsJGp)=n2;
		    valJGp(ntripletsJGp)=valJGp(ntripletsJGp)-scl.K*ajvolS_l1(i)*dJp_pj;
                case 2
                    J=mJ(l1)+epsilon_mt*dtE1;
                    dJ_v = dt*omega_p*omega_p/(dt*gamma_p+2)/linkL(l1)+ epsilon_mt/(0.5*dt*linkL(l1));
                    dJ_H = -(dt*omega_p*omega_p/(dt*gamma_p+2)+epsilon_mt/(0.5*dt));
                case 3
                    dJ_v = epsilon_in/(0.5*dt*linkL(l1));
                    dJ_H = -epsilon_in/(0.5*dt);
                    J = epsilon_in*dtE1+Js(l1);
                case 4	%all QM nodes(include the outside/inside boundary)
                    alpha_n = mun/linkL(l1); alpha_p = mup/linkL(l1);
                    beta = V(n2)-V(n1)+H(l1)*linkL(l1);
                    J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                        +Jc(alpha_p,beta,'p',p(n1),p(n2));
                    J = J_diff + epsilon_qm*dtE1;
                    dJ_v = dJndV(alpha_n,beta,n(n1),n(n2))...
                          +dJpdV(alpha_p,beta,p(n1),p(n2))+epsilon_qm/(0.5*dt*linkL(l1));
                    dJn_ni = dJdci(alpha_n,beta,'n');
                    dJn_nj = dJdcj(alpha_n,beta,'n');
                    dJp_pi = dJdci(alpha_p,beta,'p');
                    dJp_pj = dJdcj(alpha_p,beta,'p');
                    dJ_H = - dJ_v*linkL(l1);
        	    ntripletsJGn=ntripletsJGn+1;
              	    rowJGn(ntripletsJGn)=l1;
		    colJGn(ntripletsJGn)=n1;
		    valJGn(ntripletsJGn)=valJGn(ntripletsJGn)-scl.K*ajvolS_l1(i)*dJn_ni;
		    ntripletsJGn=ntripletsJGn+1;
		    rowJGn(ntripletsJGn)=l1;
		    colJGn(ntripletsJGn)=n2;
		    valJGn(ntripletsJGn)=valJGn(ntripletsJGn)-scl.K*ajvolS_l1(i)*dJn_nj;

        	    ntripletsJGp=ntripletsJGp+1;
              	    rowJGp(ntripletsJGp)=l1;
		    colJGp(ntripletsJGp)=n1;
		    valJGp(ntripletsJGp)=valJGp(ntripletsJGp)-scl.K*ajvolS_l1(i)*dJp_pi;
		    ntripletsJGp=ntripletsJGp+1;
		    rowJGp(ntripletsJGp)=l1;
		    colJGp(ntripletsJGp)=n2;
		    valJGp(ntripletsJGp)=valJGp(ntripletsJGp)-scl.K*ajvolS_l1(i)*dJp_pj;
                otherwise
                    error('undefined material');
            end
            ntripletsJGv=ntripletsJGv+1;
            rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n1;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)-scl.K*ajvolS_l1(i)*dJ_v;
	    ntripletsJGv=ntripletsJGv+1;
	    rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n2;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+scl.K*ajvolS_l1(i)*dJ_v;

            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)-scl.K*ajvolS_l1(i)*dJ_H;
            rhs_G(l1) = rhs_G(l1)-scl.K*ajvolS_l1(i)*J;
	  end
        else %mix  outside/inside QM links interface, inside QM links
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,n3+3)+epsilon_qm*dtE1);
            ntripletsJGv=ntripletsJGv+1;
            rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n1;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)- scl.K*epsilon_qm*linkS(l1)/(0.5*dt*linkL(l1)); 
	    ntripletsJGv=ntripletsJGv+1;
	    rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n2;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+ scl.K*epsilon_qm*linkS(l1)/(0.5*dt*linkL(l1));

            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+ scl.K*epsilon_qm*linkS(l1)/(0.5*dt);
        end
     end
    
JG_v=sparse(rowJGv(1:ntripletsJGv),colJGv(1:ntripletsJGv),valJGv(1:ntripletsJGv),Nlink,Nnode);
JG_H=sparse(rowJGH(1:ntripletsJGH),colJGH(1:ntripletsJGH),valJGH(1:ntripletsJGH),Nlink,Nlink);
if ~isempty(eqnSemiNodes)
JG_n=sparse(rowJGn(1:ntripletsJGn),colJGn(1:ntripletsJGn),valJGn(1:ntripletsJGn),Nlink,Nnode);
JG_p=sparse(rowJGp(1:ntripletsJGp),colJGp(1:ntripletsJGp),valJGp(1:ntripletsJGp),Nlink,Nnode);
else
JG_n=sparse(Nlink,Nnode);JG_p=sparse(Nlink,Nnode);
end
clear rowJGv colJGv valJGv rowJGH colJGH valJGH rowJGn colJGn valJGn rowJGp colJGp valJGp;
JG_v = JG_v(eqnLinks,eqnNodes);
JG_n = JG_n(eqnLinks,eqnSemiNodes);
JG_p = JG_p(eqnLinks,eqnSemiNodes);
JG_H  = JG_H(eqnLinks,eqnLinks);
JG = [JG_v,JG_n,JG_p,JG_H];
%save('dump/JG.mat','JG_v','JG_n','JG_p','JG_H');
%error('Just stop here');
clear JG_v JG_n JG_p JG_H;
    
    display(['time for matrix collection,G matrix:']);
    toc;
    tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rhs_F = rhs_F(eqnNodes);
    rhs_Kn = rhs_Kn(eqnSemiNodes);
    rhs_Kp = rhs_Kp(eqnSemiNodes);
    rhs_G = rhs_G(eqnLinks);

      
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
    dX = mylusovle(Jacob,rhs,1);
    display(['time for solving linear equation:']);
    toc;
    dV = dX(1:NeqnNodes);
    dn = dX(NeqnNodes+1:NeqnNodes+NeqnSemiNodes);
    dp = dX(NeqnNodes+NeqnSemiNodes+1:NeqnNodes+2*NeqnSemiNodes);
    dH = dX(NeqnNodes+2*NeqnSemiNodes+1:NeqnNodes+2*NeqnSemiNodes+Nl);
    
    normRes = norm(rhs);
    normRes_pre = normRes;

    itNr = itNr+1;
    
    V(eqnNodes) = V(eqnNodes)+dV;
    n(eqnSemiNodes) = n(eqnSemiNodes)+dn;
    p(eqnSemiNodes) = p(eqnSemiNodes)+dp;
    H(eqnLinks) = H(eqnLinks)+dH;

    dtV  = 2 * (V - V_p)/dt - dtV1;
    dtH  = 2 * (H - H_p)/dt - dtH1;
    dtn  = 2 * (n - n_p)/dt - dtn1;
    dtp  = 2 * (p - p_p)/dt - dtp1;

    A(eqnLinks) = A_p(eqnLinks) + (H(eqnLinks) + H_p(eqnLinks))*dt/2;
    currt = 0;


    currb = currentStd([2,round(ky/2),round(kz/2)],V,n,p,H,dtV,dtH,1);

    normUpdate = max([normRes,norm(dV)/norm(V(eqnNodes)),norm(dH)/norm(H(eqnLinks)),abs(currt)/abs(currb)]);

%    normUpdate = max([norm(dV)/norm(V(eqnNodes)),norm(dA)/norm(A(eqnLinks)),norm(dB)/norm(B(eqnLinks))]);

    display(['  Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
        '; ','normUpdate:',num2str(normUpdate),'; ','normdV:',num2str(norm(dV))]);

    maxmJ=max([mJ]);
    display(['max mJ:',num2str(maxmJ)]);
    clear dX dV dn dp dH;
    
end    
toc(tStart);   

display(['  End TD simulation relaxation2c per step.']);
