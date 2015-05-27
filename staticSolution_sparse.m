function [Vs,ns,ps] = staticSolution()

%%% static solution of Poisson and drift-diffusion (DD) equations 
%%% use Newton's method
%%% only electron density n is used (p is all zero at this moment)

display('Start static solution');

global sigma epsilon_in epsilon_sd epsilon_mt;
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
global bndLinks;
global metalNodes;
global ni;

global eqnNodes eqnSemiNodes;


%%%%%%% generate initial guess (equilibrium solution) %%%%%%%%%%%%%
Vs = zeros(Nnode,1);
dV_semiContacts = deltaV(doping(dirNodes));
Vs(dirNodes) = 1*dcVolDirNodes-dV_semiContacts;
ns = zeros(Nnode,1);
ps = zeros(Nnode,1); % only ns is used, ps is set to be zero at this version

isseminodes             = false(Nnode,1);
isseminodes(semiNodes)  = true;


for iisn=1:Nnode
    if isseminodes(iisn)
       if doping(iisn) > 0
         ns(iisn) = doping(iisn);
         ps(iisn) = ni^2./doping(iisn);
       elseif doping(iisn) < 0
         ns(iisn) = -ni^2./doping(iisn);
         ps(iisn) = -doping(iisn);
       elseif doping(iisn) == 0
         ns(iisn) = ni;
         ps(iisn) = ni;
       end
    end
end


eqnNodes = setdiff((1:Nnode)',dirNodes);
eqnSemiNodes = setdiff(semiNodes,[dirSemiNodes;dirNodes]);


%%% numbers of Poisson and DD equations
Nn = length(eqnNodes);
Nn_semi = length(eqnSemiNodes);

%%%%%%%% start Newton's iteration %%%%%%%%
newtonTol = 1e-8;
updateTol = 1e-10;
maxNewtonIt = 10;
Eps = [epsilon_sd,epsilon_mt,epsilon_in];
Sgm = [0,sigma,0];
normRes = 1;
normUpdate = 1;
itNr = 0;

tStart=tic;
while itNr < maxNewtonIt && normUpdate > updateTol
    tic;
    JF_v = zeros(Nnode,Nnode); JF_n = zeros(Nnode,Nnode); 
    JKn_v = zeros(Nnode,Nnode); JKn_n = zeros(Nnode,Nnode); 
    rhs_F = zeros(Nnode,1); rhs_Kn = zeros(Nnode,1); 
%    JF_p = []; JKn_p = []; JKp_v = []; JKp_n = [];JKp_p = []; rhs_Kp = [];
    JF_p = zeros(Nnode,Nnode); JKn_p = zeros(Nnode,Nnode); JKp_v = zeros(Nnode,Nnode);
    JKp_n = zeros(Nnode,Nnode);JKp_p = zeros(Nnode,Nnode); rhs_Kp = zeros(Nnode,1);

     ntripletsJFv=48*Nnode;
     rowJFv=zeros(ntripletsJFv,1);
     colJFv=zeros(ntripletsJFv,1);
     valJFv=zeros(ntripletsJFv,1);
     ntripletsJFv=0;

     ntripletsJFn=24*Nnode;
     rowJFn=zeros(ntripletsJFn,1);
     colJFn=zeros(ntripletsJFn,1);
     valJFn=zeros(ntripletsJFn,1);
     ntripletsJFn=0;

     ntripletsJKnv=12*Nn_semi;
     rowJKnv=zeros(ntripletsJKnv,1);
     colJKnv=zeros(ntripletsJKnv,1);
     valJKnv=zeros(ntripletsJKnv,1);
     ntripletsJKnv=0;

     ntripletsJKnn=12*Nn_semi;
     rowJKnn=zeros(ntripletsJKnn,1);
     colJKnn=zeros(ntripletsJKnn,1);
     valJKnn=zeros(ntripletsJKnn,1);
     ntripletsJKnn=0;

     ntripletsJKpv=12*Nn_semi;
     rowJKpv=zeros(ntripletsJKpv,1);
     colJKpv=zeros(ntripletsJKpv,1);
     valJKpv=zeros(ntripletsJKpv,1);
     ntripletsJKpv=0;

     ntripletsJKpp=12*Nn_semi;
     rowJKpp=zeros(ntripletsJKpp,1);
     colJKpp=zeros(ntripletsJKpp,1);
     valJKpp=zeros(ntripletsJKpp,1);
     ntripletsJKpp=0;

   % Gauss's law  %integral(epsilon*laplace V-rho)=0
   % current-continuity div*J=0
    for n1 = eqnNodes.'
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        if any(ajvolM_n1 == 1)  % node attached with semiconductor 
            dV1 = deltaV(doping(n1));
            if all(ajvolM_n1 == 1) % bulk of semiconductor %%
                rhs_F(n1) = sum(epsilon_sd*linkS(ajlk_n1).*(-(Vs(ajnd_n1)-Vs(n1)))./linkL(ajlk_n1))...
                            +(ns(n1)-ps(n1)-doping(n1))*nodeV(n1);
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    ajvol_lk = linkVolumes{lk}(1,:);
                    ajvolS_lk = linkVolumes{lk}(2,:);
                    ajvolM_lk = volumeM(ajvol_lk);
                    coefV = sum(Eps(ajvolM_lk).*ajvolS_lk)/linkL(lk);
     		    	
		    ntripletsJFv=ntripletsJFv+1;
	  	    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n1;
		    valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+coefV;
	 	    ntripletsJFv=ntripletsJFv+1;
		    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n2;
		    valJFv(ntripletsJFv)=-coefV;
                end
		ntripletsJFn=ntripletsJFn+1;
	  	rowJFn(ntripletsJFn)=n1;
		colJFn(ntripletsJFn)=n1;
		valJFn(ntripletsJFn)=nodeV(n1);

		ntripletsJFp=ntripletsJFp+1;
	  	rowJFp(ntripletsJFp)=n1;
		colJFp(ntripletsJFp)=n1;
		valJFp(ntripletsJFp)= -nodeV(n1);
           
            elseif any(ajvolM_n1 == 2)  % semi-metal interface or triple points 
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
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
                                dJp_pj = dJdcj(alpha_p,beta,'p');
				ntripletsJFn=ntripletsJFn+1;
			  	rowJFn(ntripletsJFn)=n1;
				colJFn(ntripletsJFn)=n2;
				valJFn(ntripletsJFn)=valJFn(ntripletsJFn)+ttt*ajvolS_lk(j)*dJn_nj;

				ntripletsJFp=ntripletsJFp+1;
	  			rowJFp(ntripletsJFp)=n1;
				colJFp(ntripletsJFp)=n1;
				valJFp(ntripletsJFp)=valJFp(ntripletsJFp)+ttt*ajvolS_lk(j)*dJp_pj;  
	   		    case 2
                                J = sigma*(Vs(n1)-Vs(n2)+dV1-dV2)/linkL(lk);
                                dJv = sigma/linkL(lk);  
                            case 3
                                J = 0;
                                dJv = 0;
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

                        rhs_F(n1) = rhs_F(n1)+ajvolS_lk(j)*J;
                    end
                end                 
            elseif any(ajvolM_n1 == 3) % semi-insulator interface
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
		    ntripletsJFv=ntripletsJFv+1;
	  	    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n1;
		    valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+coefV;
	 	    ntripletsJFv=ntripletsJFv+1;
		    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n2;
		    valJFv(ntripletsJFv)=-coefV;
                    dV2_vec = [0,1,1/2];
                    rhs_F(n1) = rhs_F(n1)+sum(Eps(ajvolM_lk).*...
                        (-(Vs(n2)-Vs(n1)+dV2_vec(ajvolM_lk).*dV2)).*ajvolS_lk)/linkL(lk);
                end
                semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
		ntripletsJFn=ntripletsJFn+1;
	  	rowJFn(ntripletsJFn)=n1;
		colJFn(ntripletsJFn)=n1;
		valJFn(ntripletsJFn)=semiV;

		ntripletsJFp=ntripletsJFp+1;
	  	rowJFp(ntripletsJFp)=n1;
		colJFp(ntripletsJFp)=n1;
		valJFp(ntripletsJFp)= -semiV;
                rhs_F(n1) = rhs_F(n1)+(ns(n1)-ps(n1)-doping(n1))*semiV;
            else
                error('Incorrect node assignment');
            end
        else                   % node attached no semiconductor
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

                if any(ajvolM_n1 == 2)
                 tcoem = 0;
                else 
                 tcoem = 1;
                end

                coefV = sum((tcoem*Eps(ajvolM_lk)+Sgm(ajvolM_lk)).*ajvolS_lk)/linkL(lk);
		    ntripletsJFv=ntripletsJFv+1;
	  	    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n1;
		    valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+coefV;
	 	    ntripletsJFv=ntripletsJFv+1;
		    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n2;
		    valJFv(ntripletsJFv)=-coefV;
                dV2_vec = [0,1,1/2];
                rhs_F(n1) = rhs_F(n1)+sum((tcoem*Eps(ajvolM_lk)+Sgm(ajvolM_lk)).*...
                        (-(Vs(n2)-Vs(n1)+dV2_vec(ajvolM_lk).*dV2)).*ajvolS_lk)/linkL(lk);
            end
        end
       % if any(ajvolM_n1 == 2) % scaling some equations with condutivity for matrix balance
       %     JF_v(n1,:) = JF_v(n1,:)/sigma;
       %     JF_n(n1,:) = JF_n(n1,:)/sigma;
       %     JF_p(n1,:) = JF_p(n1,:)/sigma;
       %     rhs_F(n1) = rhs_F(n1)/sigma;
       % end
    end
JF_v=sparse(rowJFv(1:ntripletsJFv),colJFv(1:ntripletsJFv),valJFv(1:ntripletsJFv),Nnode,Nnode);
JF_n=sparse(rowJFn(1:ntripletsJFn),colJFn(1:ntripletsJFn),valJFn(1:ntripletsJFn),Nnode,Nnode);
JF_p=sparse(rowJFp(1:ntripletsJFp),colJFp(1:ntripletsJFp),valJFp(1:ntripletsJFp),Nnode,Nnode);
JF_v = JF_v(eqnNodes,eqnNodes);
JF_n = JF_n(eqnNodes,eqnSemiNodes);
JF_p = JF_p(eqnNodes,eqnSemiNodes);
display(['time for matrix collection,F matrix:']);
toc;
tic;

    %%%%%%%%% build Jacobian and rhs of K %%%%%%%%%%%%%!for semiconductor
    if ~isempty(eqnSemiNodes)
        for n1 = eqnSemiNodes.'
            ajlk_n1 = nodeLinks{n1}(1,:);
            ajnd_n1 = nodeLinks{n1}(2,:);
%             ajvol_n1 = nodeVolumes{n1}(1,:);
            for i = 1:length(ajlk_n1)
                n2 = ajnd_n1(i);
                lk = ajlk_n1(i);
                ajvol_lk = linkVolumes{lk}(1,:);
                ajvolS_lk = linkVolumes{lk}(2,:);
                ajvolM_lk = volumeM(ajvol_lk);
                semiS = sum(ajvolS_lk(ajvolM_lk == 1));
                alpha_n = mun/linkL(lk);  alpha_p = mup/linkL(lk);
                beta = Vs(n2)-Vs(n1);
                
                dJn_v = dJndV(alpha_n,beta,ns(n1),ns(n2));
	    	ntripletsJKnv=ntripletsJKnv+1;
		rowJKnv(ntripletsJKnv)=n1;
		colJKnv(ntripletsJKnv)=n1;
		valJKnv(ntripletsJKnv)=valJKnv(ntripletsJKnv)+semiS*dJn_v;
	    	ntripletsJKnv=ntripletsJKnv+1;
		rowJKnv(ntripletsJKnv)=n1;
		colJKnv(ntripletsJKnv)=n2;
		valJKnv(ntripletsJKnv)=valJKnv(ntripletsJKnv)-semiS*dJn_v;
                dJp_v = dJpdV(alpha_p,beta,ps(n1),ps(n2));
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
                
                Jn = Jc(alpha_n,beta,'n',ns(n1),ns(n2));
                Jp = Jc(alpha_p,beta,'p',ps(n1),ps(n2));
                rhs_Kn(n1) = rhs_Kn(n1)+Jn*semiS;
                rhs_Kp(n1) = rhs_Kp(n1)+Jp*semiS;
            end
        end
    end

if ~isempty(eqnSemiNodes)
JKn_v=sparse(rowJKnv(1:ntripletsJKnv),colJKnv(1:ntripletsJKnv),valJKnv(1:ntripletsJKnv),Nnode,Nnode);
JKn_n=sparse(rowJKnn(1:ntripletsJKnn),colJKnn(1:ntripletsJKnn),valJKnn(1:ntripletsJKnn),Nnode,Nnode);
JKp_v=sparse(rowJKpv(1:ntripletsJKpv),colJKpv(1:ntripletsJKpv),valJKpv(1:ntripletsJKpv),Nnode,Nnode);
JKp_p=sparse(rowJKpp(1:ntripletsJKpp),colJKpp(1:ntripletsJKpp),valJKpp(1:ntripletsJKpp),Nnode,Nnode);
JKn_v = JKn_v(eqnSemiNodes,eqnNodes);
JKn_n = JKn_n(eqnSemiNodes,eqnSemiNodes);
JKn_p = [];
JKp_v = JKp_v(eqnSemiNodes,eqnNodes);
JKp_p = JKp_p(eqnSemiNodes,eqnSemiNodes);
JKp_n = [];
else
JKn_v=[];JKn_H=[];JKn_p=[];JKn_n=[];
JKp_v=[];JKp_H=[];JKp_n=[];JKp_p=[];
end
display(['time for matrix collection,K matrix:']);
toc;
tic;

rhs_F = rhs_F(eqnNodes);
rhs_Kn =rhs_Kn(eqnSemiNodes);
rhs_Kp =rhs_Kp(eqnSemiNodes);
Jacob = [JF_v,JF_n,JF_p;JKn_v,JKn_n,JKn_p;JKp_v,JKp_n,JKp_p];
rhs = -[rhs_F;rhs_Kn;rhs_Kp];

    clear JF_v JF_n JF_p;
    clear JKn_v JKn_n JKn_p;
    clear JKp_v JKp_n JKp_p;
    clear rhs_F rhs_Kn rhs_Kp;

%    invJacob = Jacob^(-1);
%    dX = invJacob*rhs;
    dX = mylusovle(Jacob,rhs,1);
    display(['time for solving linear equation:']);
    toc;
    dV = dX(1:Nn);
    dn = dX(Nn+1:Nn+Nn_semi);
    dp = dX(Nn+Nn_semi+1:Nn+2*Nn_semi);
    normRes = norm(rhs);
%    normUpdate = max(norm(dV)/norm(Vs(eqnNodes)));
    itNr = itNr+1;
    
    Vs(eqnNodes) = Vs(eqnNodes)+dV;
    ns(eqnSemiNodes) = ns(eqnSemiNodes)+dn;
    ps(eqnSemiNodes) = ps(eqnSemiNodes)+dp;    

    currb=currentS([2,round(ky/2),round(kz/2)],Vs,ns,ps,1);

    currt=currentS([1,round(ky/2),round(kz/2)],Vs,ns,ps,1)-currentS([kx,round(ky/2),round(kz/2)],Vs,ns,ps,1)+currentS([round(kx/2),1,round(kz/2)],Vs,ns,ps,2)-currentS([round(kx/2),ky,round(kz/2)],Vs,ns,ps,2)+currentS([round(kx/2),round(ky/2),1],Vs,ns,ps,3)-currentS([round(kx/2),round(ky/2),kz],Vs,ns,ps,3);
	
%	normUpdate = max(norm(dV)/norm(Vs(eqnNodes)),norm(currb-curra)/norm(curra));

    normUpdate = max([normRes,norm(dV)/norm(Vs(eqnNodes)),abs(currt)/abs(currb)]);

	display(['Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
	        '; ','normUpdate:',num2str(normUpdate),'; ','normdV:',num2str(norm(dV))]);

	clear dX dV dn dp;
    	clear Jacob rhs;
end  
% 
% Vs = Vs.*scl.Vt;
% ns = ns.*scl.ni;
% ps = ps.*scl.ni;

toc(tStart);   

display('End static solution');
