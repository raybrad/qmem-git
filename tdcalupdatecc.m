function [dtV,dtn,dtp,dtH,itNr] = tdcalupdatecc(V,n,p,A,H,mJ_p,Js,Yp,Zp)

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

global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global isqmnodes;

global currdlink;

global ehcrf;

global isqmvolm;


%%%%%%% initial guess  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y    = Yp;
dtn  = zeros(Nnode,1); 
dtp  = zeros(Nnode,1);
Z    = Zp;


nobndNodes = setdiff((1:Nnode)',bndNodes);

eqnNodes  = setdiff(bndNodes,[dirNodes]);
eqnNodes1 = setdiff(nobndNodes,[dirNodes]); 
eqnSemiNodes = setdiff(semiNodes,[dirSemiNodes;dirNodes;sQMnodes]);

eqnLinks = setdiff((1:Nlink)',[bndLinks]);

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


%%%% Lorezentz gauge%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par               par     1                      %
%---- V= div*A  --> ----V+ ---- integral0(A)dS=0   %
%par t              par t   Vol                    %
%%%%%%% initial guess  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n1 = eqnNodes1.'
    ajlk_n1 = nodeLinks{n1}(1,:);
    ajnd_n1 = nodeLinks{n1}(2,:);
    ajvol_n1 = nodeVolumes{n1}(1,:);
    ajvolV_n1 = nodeVolumes{n1}(2,:);
    ajvolM_n1 = volumeM(ajvol_n1);
    for qmi = 1:length(ajvol_n1)
        if isqmvolm(ajvol_n1(qmi)) 
           ajvolM_n1(qmi)=4;
        end
    end
    sign_n1 = sign(ajnd_n1-n1);
    coef_n1 = linkS(ajlk_n1)./nodeV(n1);
    Y(n1)  = -sign_n1.*coef_n1'*A(ajlk_n1);
    epsnode = sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1);
    Y(n1) = Y(n1) / (scl.K * epsnode);
end

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
           if isSQMlinks(lk) == 0
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
           elseif isSQMlinks(lk) == 1
             dtn(n1)   = dtn(n1) + ehcrf     * linkS(lk)*currdlink(lk,4);
             dtp(n1)   = dtp(n1) - (1-ehcrf) * linkS(lk)*currdlink(lk,4);
           elseif isSQMlinks(lk) == 2
             dtn(n1)   = dtn(n1) - ehcrf     * linkS(lk)*currdlink(lk,4);
             dtp(n1)   = dtp(n1) + (1-ehcrf) * linkS(lk)*currdlink(lk,4);
           elseif isSQMlinks(lk) == 3
             dtn(n1)   = dtn(n1) + ehcrf     * linkS(lk)*currdlink(lk,5);
             dtp(n1)   = dtp(n1) - (1-ehcrf) * linkS(lk)*currdlink(lk,5);
           elseif isSQMlinks(lk) == 4
             dtn(n1)   = dtn(n1) - ehcrf     * linkS(lk)*currdlink(lk,5);
             dtp(n1)   = dtp(n1) + (1-ehcrf) * linkS(lk)*currdlink(lk,5);
           elseif isSQMlinks(lk) == 5
             dtn(n1)   = dtn(n1) + ehcrf     * linkS(lk)*currdlink(lk,6);
             dtp(n1)   = dtp(n1) - (1-ehcrf) * linkS(lk)*currdlink(lk,6);
           elseif isSQMlinks(lk) == 6
             dtn(n1)   = dtn(n1) - ehcrf     * linkS(lk)*currdlink(lk,6);
             dtp(n1)   = dtp(n1) + (1-ehcrf) * linkS(lk)*currdlink(lk,6);
           elseif isSQMlinks(lk) > 6
             dtn(n1)   = dtn(n1) + 0;
             dtp(n1)   = dtp(n1) - 0;
           end
       end
       semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
       if semiV ~=0
          dtn(n1) = dtn(n1)/semiV;
          dtp(n1) = dtp(n1)/semiV;
       end
   end
end
    
while itNr < maxNewtonIt && normUpdate > updateTol

     rhs_F = zeros(Nnode,1); 
     rhs_G = zeros(Nlink,1);

     ntripletsJFv=12*Nnode;
     rowJFv=zeros(ntripletsJFv,1);
     colJFv=zeros(ntripletsJFv,1);
     valJFv=zeros(ntripletsJFv,1);
     ntripletsJFv=0;

     ntripletsJFH=6*Nnode;
     rowJFH=zeros(ntripletsJFH,1);
     colJFH=zeros(ntripletsJFH,1);
     valJFH=zeros(ntripletsJFH,1);
     ntripletsJFH=0;

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
     
     for k=1:NeqnNodes
        n1 = eqnNodes(k);
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        sign_n1 = sign(ajnd_n1-n1);
        if isDirSemiNodes(n1), dV1 = deltaV(doping(n1)); else dV1 = 0; end
        for i = 1:length(ajlk_n1)
            n2 = ajnd_n1(i);
            lk = ajlk_n1(i);
            ajvol_lk = linkVolumes{lk}(1,:);
            ajvolS_lk = linkVolumes{lk}(2,:);
            ajvolM_lk = volumeM(ajvol_lk);
            for qmi = 1:length(ajvol_lk)
             if isqmvolm(ajvol_lk(qmi))
               ajvolM_lk(qmi)=4;
             end
            end

            if isDirSemiNodes(n2), dV2 = deltaV(doping(n2)); else dV2 = 0; end

            dtE1 = -(Y(n2)-Y(n1))/linkL(lk)-sign_n1(i)*Z(lk);
            for j = 1:length(ajvol_lk)
               switch ajvolM_lk(j)
                 case 1
                  alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                  beta = V(n2)-V(n1) + sign_n1(i)*H(lk)*linkL(lk);
                  J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                          +Jc(alpha_p,beta,'p',p(n1),p(n2));
                  Jd   = epsilon_sd * dtE1;
                  J    = J_diff +Jd;
                  dJv  = epsilon_sd/linkL(lk);
                  dJH  =-sign_n1(i)*epsilon_sd;
                 case 2
                  J=  sign_n1(i)*mJ_p(lk) + epsilon_mt*dtE1;
                  dJv  = epsilon_mt/linkL(lk);
                  dJH  =-sign_n1(i)*epsilon_mt;
                 case 3
                  J   = epsilon_in*dtE1;
                  dJv  = epsilon_in/linkL(lk);
                  dJH  =-sign_n1(i)*epsilon_in;
                 case 4
                  alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
                  beta = V(n2)-V(n1) + sign_n1(i)*H(lk)*linkL(lk);
                  J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                          +Jc(alpha_p,beta,'p',p(n1),p(n2));
                  Jd   = epsilon_qm * dtE1;
                  J    = J_diff +Jd;
                  dJv  = epsilon_qm/linkL(lk);
                  dJH  =-sign_n1(i)*epsilon_qm;
                 otherwise
                  error('undefined material');
               end
                rhs_F(n1) = rhs_F(n1)+ajvolS_lk(j)*J;

		ntripletsJFv=ntripletsJFv+1;
		rowJFv(ntripletsJFv)=n1;
		colJFv(ntripletsJFv)=n1;
		valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+ajvolS_lk(j)*dJv;


               if ~isBndNodes(n2) 
		 ntripletsJFv=ntripletsJFv+1;
		 rowJFv(ntripletsJFv)=n1;
		 colJFv(ntripletsJFv)=n2;
		 valJFv(ntripletsJFv)=valJFv(ntripletsJFv)-0;
               else
	 	 ntripletsJFv=ntripletsJFv+1;
		 rowJFv(ntripletsJFv)=n1;
		 colJFv(ntripletsJFv)=n2;
		 valJFv(ntripletsJFv)=valJFv(ntripletsJFv)-ajvolS_lk(j)*dJv;
               end

	 	 ntripletsJFH=ntripletsJFH+1;
		 rowJFH(ntripletsJFH)=n1;
		 colJFH(ntripletsJFH)=lk;
		 valJFH(ntripletsJFH)=valJFH(ntripletsJFH)+ajvolS_lk(j)*dJH;
            end
        end
     end

    %%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
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
        CC = zeros(1,Nlink);

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

        %%%%%% source current %%%%%%%%%%%
        dtE1 = -(Y(n2)-Y(n1))/linkL(l1) - Z(l1);
        if isSQMlinks(l1) == 0
          E1 = -(V(n2)-V(n1))/linkL(l1)-H(l1);
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
                    dJ_v = epsilon_sd/linkL(l1);
                    dJ_H = -epsilon_sd; 
                case 2
                    dJ_v = epsilon_mt/linkL(l1);
                    dJ_H = -epsilon_mt;
                    J =  mJ_p(l1) + epsilon_mt*dtE1;
                case 3
                    dJ_v = epsilon_in/linkL(l1);
                    dJ_H = -epsilon_in;
                    J = epsilon_in*dtE1+Js(l1);
                case 4
                    alpha_n = mun/linkL(l1); alpha_p = mup/linkL(l1);
                    beta = V(n2)-V(n1)+H(l1)*linkL(l1);
                    J_diff = Jc(alpha_n,beta,'n',n(n1),n(n2))...
                        +Jc(alpha_p,beta,'p',p(n1),p(n2));
                    J = J_diff + epsilon_qm*dtE1;
                    dJ_v = epsilon_qm/linkL(l1);
                    dJ_H = -epsilon_qm;
                otherwise
                    error('undefined material');
            end

            if ~isBndNodes(n1)
              ntripletsJGv=ntripletsJGv+1;
              rowJGv(ntripletsJGv)=l1;
	      colJGv(ntripletsJGv)=n1;
	      valJGv(ntripletsJGv)=valJGv(ntripletsJGv)-0;
            else
              ntripletsJGv=ntripletsJGv+1;
              rowJGv(ntripletsJGv)=l1;
	      colJGv(ntripletsJGv)=n1;
	      valJGv(ntripletsJGv)=valJGv(ntripletsJGv)-scl.K*ajvolS_l1(i)*dJ_v;
            end

            if ~isBndNodes(n2)
	      ntripletsJGv=ntripletsJGv+1;
	      rowJGv(ntripletsJGv)=l1;
	      colJGv(ntripletsJGv)=n2;
	      valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+0;
            else 
	      ntripletsJGv=ntripletsJGv+1;
	      rowJGv(ntripletsJGv)=l1;
	      colJGv(ntripletsJGv)=n2;
	      valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+scl.K*ajvolS_l1(i)*dJ_v;
            end

            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)-scl.K*ajvolS_l1(i)*dJ_H;
            rhs_G(l1) = rhs_G(l1)-scl.K*ajvolS_l1(i)*J;
          end
        elseif isSQMlinks(l1) == 1
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,4)+epsilon_qm*dtE1);
            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+scl.K*epsilon_qm*linkS(l1);
        elseif isSQMlinks(l1) == 2
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,4)+epsilon_qm*dtE1);
            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+scl.K*epsilon_qm*linkS(l1);
        elseif isSQMlinks(l1) == 3
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,5)+epsilon_qm*dtE1);
            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+scl.K*epsilon_qm*linkS(l1);
        elseif isSQMlinks(l1) == 4
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,5)+epsilon_qm*dtE1);
            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+scl.K*epsilon_qm*linkS(l1);
        elseif isSQMlinks(l1) == 5
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,6)+epsilon_qm*dtE1);
            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+scl.K*epsilon_qm*linkS(l1);
        elseif isSQMlinks(l1) == 6
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,6)+epsilon_qm*dtE1);
            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+scl.K*epsilon_qm*linkS(l1);
        elseif isSQMlinks(l1) == 7
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,4)+epsilon_qm*dtE1);
            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+scl.K*epsilon_qm*linkS(l1);
        elseif isSQMlinks(l1) == 8
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,5)+epsilon_qm*dtE1);
            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+scl.K*epsilon_qm*linkS(l1);
        elseif isSQMlinks(l1) == 9
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,6)+epsilon_qm*dtE1);
            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+scl.K*epsilon_qm*linkS(l1);
        end
     end
    
    JF_v=sparse(rowJFv(1:ntripletsJFv),colJFv(1:ntripletsJFv),valJFv(1:ntripletsJFv),Nnode,Nnode);%ntripletsJFv,ntripletsJFv);
    JF_H=sparse(rowJFH(1:ntripletsJFH),colJFH(1:ntripletsJFH),valJFH(1:ntripletsJFH),Nnode,Nlink);%ntripletsJFH,ntripletsJFH);
    clear rowJFv colJFv valJFv rowJFH colJFH valJFH;
    JF_v = JF_v(eqnNodes,eqnNodes);
    JF_H  = JF_H(eqnNodes,eqnLinks);
    
    JG_v=sparse(rowJGv(1:ntripletsJGv),colJGv(1:ntripletsJGv),valJGv(1:ntripletsJGv),Nlink,Nnode);%ntripletsJGv,ntripletsJGv);
    JG_H=sparse(rowJGH(1:ntripletsJGH),colJGH(1:ntripletsJGH),valJGH(1:ntripletsJGH),Nlink,Nlink);%ntripletsJGH,ntripletsJGH);
    clear rowJGv colJGv valJGv rowJGH colJGH valJGH;
    JG_v = JG_v(eqnLinks,eqnNodes);
    JG_H  = JG_H(eqnLinks,eqnLinks);
    
    rhs_F = rhs_F(eqnNodes);
    rhs_G = rhs_G(eqnLinks);

    JF = [JF_v,JF_H];
    JG = [JG_v,JG_H];

    clear JF_v JF_H JG_v JG_H;

    Jacob = [JF;JG];
    rhs = -[rhs_F;rhs_G];
  
    clear JF JG rhs_F rhs_G;
    dX = Jacob\rhs;
    dY = dX(1:NeqnNodes);
    dZ = dX(NeqnNodes+1:NeqnNodes+Nl);
    
    normRes = norm(rhs);
    normRes_pre = normRes;

    itNr = itNr+1;
    
    Y(eqnNodes) = Y(eqnNodes)+dY;
    Z(eqnLinks) = Z(eqnLinks)+dZ;

    currt=currentStd([1,round(ky/2),round(kz/2)],V,n,p,H,Y,Z,1)-currentStd([kx,round(ky/2),round(kz/2)],V,n,p,H,Y,Z,1)+currentStd([round(kx/2),1,round(kz/2)],V,n,p,H,Y,Z,2)-currentStd([round(kx/2),ky,round(kz/2)],V,n,p,H,Y,Z,2)+currentStd([round(kx/2),round(ky/2),1],V,n,p,H,Y,Z,3)-currentStd([round(kx/2),round(ky/2),kz],V,n,p,H,Y,Z,3);

    currb = currentStd([2,round(ky/2),round(kz/2)],V,n,p,H,Y,Z,1);

    normUpdate = max([normRes,norm(dY)/norm(Y(eqnNodes)),norm(dZ)/norm(Z(eqnLinks)),abs(currt)/abs(currb)]);
    
    clear dX dY dZ;    
end    

dtV = Y;
dtH = Z;
display(['  Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),'; ','normUpdate:',num2str(normUpdate)]);
