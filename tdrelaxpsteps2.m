function [V,A,H,mJ_0,mJ_1,mJ_2,mJ_1p,mJ_2p,Efield_p,dtV,dtH] = tdrelaxpstep2(V_p,A_p,H_p,Js,mJ_0p,mJ_1p,mJ_2p,mJ_1pp,mJ_2pp,Efield_pp,dtV,dtH,dt)

display(['  Start TD simulation relaxation2 per step:']);

global sigma epsilon_in;
global epsilon_mt;
%
global omega_p gamma_p;
global lepsr_1 lomega_1 lgamma_1 lepsr_2 lomega_2 lgamma_2;
global scl;
global nodes links contacts;
global Nnode Nlink;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes dirNodes;
global isBndNodes isDirNodes dcVolDirNodes acVolDirNodes;
global kx ky kz;
global bndLinks;
global metalNodes metalLinks;
global savefile;
global light_speed XZsurfLinks XYsurfLinks YZsurfLinks;
global lightdirection EsurfLinks BsurfLinks;

%Method 2
%Solving Gauss' law.
%Solving Current-continuity equations.
%Solving Maxwell-Ampere equations(already insert Lorentz gauge)
%Respecting the gauge condition is a side product of solving above set.

%%%%%%% initial guess  %%%%%%%%%%%%%
V = V_p + dtV *dt;
H = H_p + dtH *dt;
A = A_p + H_p *dt + dtH * dt^2/2;

mJ_0=mJ_0p;
mJ_1=mJ_1p;
mJ_2=mJ_2p;
Efield=Efield_pp;
Efield_p=Efield_pp;

dtV1 = dtV;
dtH1 = dtH;

eqnNodes = setdiff((1:Nnode)',dirNodes);%eqnNodes (nodes except from contact ,different from tdcalupdatec.m)   
					%dirNodes (contact part update by dtVp)

eqnLinks=setdiff((1:Nlink)',EsurfLinks);	%for E boundary condition
isBsurfLinks=false(Nlink,1);
isBsurfLinks(BsurfLinks)=true;		        %for B boundary condition
%%%%%%%% start Newton's iteration %%%%%%%%
updateTol = 1e-8;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 50;
maxLinIt = 200;
Eps = [epsilon_mt,epsilon_in];
Sgm = [sigma,0];

NeqnNodes = length(eqnNodes);
Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
normRes_pre = 1e10;
itNr = 0;
tStart=tic;
while itNr < maxNewtonIt && normUpdate > updateTol
tic;
     % construct the sparse matrix by (Row,Column,Value) 
     rhs_F = zeros(Nnode,1); 
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

    %update mJ,for plasmonic metal
    for l1 = metalLinks'
    n1 = links(l1,1);
    n2 = links(l1,2);
    %drude mJ_(n+1)=(2-gamma_p*dt)/(2+gamma_p*dt)*mJ_n+ omega_p^2*dt/(2+gamma_p*dt)*(E_n+1 + E_n)
    Efield(l1)=   -(V(n2)-V(n1))/linkL(l1)-H(l1);
    Efield_p(l1)= -(V_p(n2)-V_p(n1))/linkL(l1)-H_p(l1);
    mJ_0(l1) =(2-gamma_p*dt)/(2+gamma_p*dt)*mJ_0p(l1)+dt*omega_p*omega_p/(2+gamma_p*dt)*(Efield(l1)+Efield_p(l1));
   %alpha=(2-lomega*lomega*dt*dt)/(1+lgamma*dt);beta=(lgamma*dt-1)/(1+lgamma*dt);gamma=epsr*lomega*lomega*dt*dt/(1+lgamma*dt)
   %Lorentz mJ_(n+1)=alpha*mJ_n+beta*mJ_(n-1) + gamma * (E_n+1 -E_n-1)/2/dt 
   %or to save the storage, Lorentz mJ_(n+1)=alpha*mJ_n+beta*mJ_(n-1) + gamma* (E_n+1 -E_n)/dt 
    mJ_1(l1)=(2-lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)*mJ_1p(l1)+(lgamma_1*dt-1)/(1+lgamma_1*dt)*mJ_1pp(l1)+...
    	     (lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt)*(Efield(l1)-Efield_pp(l1));

    mJ_2(l1)=(2-lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)*mJ_2p(l1)+(lgamma_2*dt-1)/(1+lgamma_2*dt)*mJ_2pp(l1)+...
    	     (lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt)*(Efield(l1)-Efield_pp(l1));
    end
display(['time for updating metalLinks:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rhsF
   for k=1:NeqnNodes
        n1 = eqnNodes(k);
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        sign_n1 = sign(ajnd_n1-n1);

        if any(ajvolM_n1 == 1) 
            for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    ajvol_lk = linkVolumes{lk}(1,:);
                    ajvolS_lk = linkVolumes{lk}(2,:);
                    ajvolM_lk = volumeM(ajvol_lk);
 
                      dV2 = 0;
                      ttt = 1; 

                    dtE1 = -(dtV(n2)-dtV(n1))/linkL(lk)-sign_n1(i)*dtH(lk);
				
                    for j = 1:length(ajvol_lk)
                        switch ajvolM_lk(j)
                            case 1
                  		J=sign_n1(i)*(mJ_0(lk)+mJ_1(lk)+mJ_2(lk))+epsilon_mt*dtE1;
                                dJv =(dt*omega_p*omega_p/(dt*gamma_p+2)/linkL(lk)+...
				    (lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt)/linkL(lk) +...
				    (lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt)/linkL(lk)) +...
				     epsilon_mt/(0.5*dt*linkL(lk));
                                dJH = -sign_n1(i)*(dt*omega_p*omega_p/(dt*gamma_p+2) +...
				                  (lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt) +...
				                  (lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt) +...
						  epsilon_mt/(0.5*dt));
                            case 2
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
            end
        else
                for i = 1:length(ajlk_n1)

                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                        dV2 = 0;
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

                    dV2_vec = [1,1/2];
                    coefH = -coefV*linkL(lk);

		    ntripletsJFH=ntripletsJFH+1;
		    rowJFH(ntripletsJFH)=n1;
		    colJFH(ntripletsJFH)=lk;
		    valJFH(ntripletsJFH)=sign_n1(i)*coefH;

                    rhs_F(n1) = rhs_F(n1)+sum(Eps(ajvolM_lk).*(-(V(n2)-V(n1)+dV2_vec(ajvolM_lk).*dV2)/linkL(lk)-sign_n1(i)* H(lk)).*ajvolS_lk);
                end
        end
     end
JF_v=sparse(rowJFv(1:ntripletsJFv),colJFv(1:ntripletsJFv),valJFv(1:ntripletsJFv),Nnode,Nnode);
JF_H=sparse(rowJFH(1:ntripletsJFH),colJFH(1:ntripletsJFH),valJFH(1:ntripletsJFH),Nnode,Nlink);
clear rowJFv colJFv valJFv rowJFH colJFH valJFH;
JF_v = JF_v(eqnNodes,eqnNodes);
JF_H = JF_H(eqnNodes,eqnLinks);
JF = [JF_v,JF_H];
clear JF_v JF_H;
display(['time for matrix collection,F matrix:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rhsG

    %%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
    % G = cur curl A - {grad div A+grad V}+{-J+par E/par t} =0 (for link A,PI)
    %  already contains lorentz gauge
    %PI_n+1 =PI_n - G/{par G/parPI}  to get the PI
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

    for k=1:Nl 
        l1 = eqnLinks(k);
        n1 = links(l1,1);	%sta node
        n2 = links(l1,2);	%end node
        ajvol_l1 = linkVolumes{l1}(1,:);%volumeid
        ajvolS_l1 = linkVolumes{l1}(2,:);%associate surf area,not dual area, but only its own 1/4
        ajvolM_l1 = volumeM(ajvol_l1);
        CC = zeros(1,Nlink);
        GD = zeros(1,Nlink);
        bndJG_H = zeros(1,Nlink);
	ajlk_n1=[];
	ajlk_n2=[];
        %%%%%%%%%% curl curl of A %%%%%%%%%%
        ajsf_l1 = linkSurfs{l1}; %the surfLink l1 belong to(at most 4),which is in the form of {(surf index,link);(surf index,link);...}
        for i = 1:size(ajsf_l1,2)	%num of cell ,i.e. num of surfs the link is on
            ajlk = ajsf_l1(2:4,i);	% the other 3 link on that surface i
            S = linkL(l1)*linkL(ajlk(1)); %surface area the link l1 is on, linkL ~link length  
            dL = dlinkL(ajsf_l1(1,i));    %1~surf index,dLinkL~dual length of links, (Nsurf,1)
            ajlkL = linkL(ajlk)';	 % length of the other 3 link
            if n2 == links(ajlk(1),1)	%the end node of link l1 is the start node of the first surround link 
                ajlkSign = [1,-1,-1];
            elseif n2 == links(ajlk(1),2)
                ajlkSign = [-1,-1,1];
            else
                error('Incorrect link arrangement');
            end
            CC(ajlk) = dL/S*(ajlkSign.*ajlkL); %CC(1:3) other 3 link ,those link could be on the boundary,just l1 is not bnd
            CC(l1) = CC(l1)+dL/S*linkL(l1);	% link l1
            rhs_G(l1) = rhs_G(l1)+CC(ajlk)*A(ajlk);  %A(ajlk) contain info from boundary
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

	   %if (size(ajsf_l1,2)==3) %xy plane : Ay(lk)~Bx~par Ay(lk)/par z =1/c par Ay(lk)/par t no matter top or bottom
	   %    if(isBsurfLinks(l1))
	   %    rhs_G(l1)=rhs_G(l1)+0.0;
	   %    bndJG_H(l1)=0.0;
	   %    ntripletsJGH=ntripletsJGH+1;
	   %    rowJGH(ntripletsJGH)=l1;
	   %    colJGH(ntripletsJGH)=l1;
	   %    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
       	   %else
       	   %    rhs_G(l1)=rhs_G(l1)+dL*(1/light_speed*H(l1));
       	   %    bndJG_H(l1)=dL/light_speed;
	   %    ntripletsJGH=ntripletsJGH+1;
	   %    rowJGH(ntripletsJGH)=l1;
	   %    colJGH(ntripletsJGH)=l1;
	   %    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
	   %end
	   %elseif (size(ajsf_l1,2)==2) 
	   %    if(isBsurfLinks(l1))
	   %    rhs_G(l1)=rhs_G(l1)+dL*(1/light_speed*H(l1));
	   %    bndJG_H(l1)=dL/light_speed;
	   %    ntripletsJGH=ntripletsJGH+1;
	   %    rowJGH(ntripletsJGH)=l1;
	   %    colJGH(ntripletsJGH)=l1;
	   %    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
	   %    else
	   %    rhs_G(l1)=rhs_G(l1)+2*dL*(1/light_speed*H(l1));
	   %    bndJG_H(l1)=2*dL/light_speed;
	   %    ntripletsJGH=ntripletsJGH+1;
	   %    rowJGH(ntripletsJGH)=l1;
	   %    colJGH(ntripletsJGH)=l1;
	   %    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
	   %    end
	   %end
	       rhs_G(l1)=rhs_G(l1)+(1-isBsurfLinks(l1) + 3 - size(ajsf_l1,2))* dL*(1/light_speed*H(l1));
	       bndJG_H(l1)=(1-isBsurfLinks(l1) + 3 - size(ajsf_l1,2))*dL/light_speed;
	       ntripletsJGH=ntripletsJGH+1;
	       rowJGH(ntripletsJGH)=l1;
	       colJGH(ntripletsJGH)=l1;
	       valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
        %%%%%%%% gradient divergence of A and gradient of V %%%%%%%%
	%innerlinks has both not BndNodes(n1) and not BndNodes(n2)
	%links connect to bndNodes,only has half
	%if both nodes are on surface,no Lorentz
        if ~isBndNodes(n1)
            ajlk_n1 = nodeLinks{n1}(1,:); % sta node connected to at most 6 links
            ajnd_n1 = nodeLinks{n1}(2,:); % sta node connected to at most 6 nodes
            sign_n1 = sign(ajnd_n1-n1);
            coef_n1 = linkS(ajlk_n1)./nodeV(n1)*linkS(l1)/linkL(l1); %linkS(Nlink,1)~ dual surface area,nodeV(Nnode,1)
            GD(ajlk_n1) = GD(ajlk_n1)-sign_n1.*coef_n1';
            ajvol_n1 = nodeVolumes{n1}(1,:);
            ajvolV_n1 = nodeVolumes{n1}(2,:);
            ajvolM_n1 = volumeM(ajvol_n1);
            coefV_n1 = scl.K*(sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1))*linkS(l1)/linkL(l1);
	
            ntripletsJGv=ntripletsJGv+1;
	    rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n1;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+coefV_n1/(0.5*dt);

            rhs_G(l1) = rhs_G(l1)-(-sign_n1.*coef_n1')*A(ajlk_n1); %gradient divergence A,if it's connected to a bndNode,linkS
	    							   %could be on bndsurf
            rhs_G(l1) = rhs_G(l1)+coefV_n1*dtV(n1);		   %dV/dt
        end

        if ~isBndNodes(n2)
            ajlk_n2 = nodeLinks{n2}(1,:);%6 neighbor links around node n2
            ajnd_n2 = nodeLinks{n2}(2,:);%6 neighboring nodes around node n2 
            sign_n2 = sign(ajnd_n2-n2);
            coef_n2 = linkS(ajlk_n2)./nodeV(n2)*linkS(l1)/linkL(l1);
            GD(ajlk_n2) = GD(ajlk_n2)+sign_n2.*coef_n2';
            ajvol_n2 = nodeVolumes{n2}(1,:);	%volumeid
            ajvolV_n2 = nodeVolumes{n2}(2,:);	%associate 1/8 volume
            ajvolM_n2 = volumeM(ajvol_n2);
            coefV_n2 = scl.K*(sum(Eps(ajvolM_n2).*ajvolV_n2)/nodeV(n2))*linkS(l1)/linkL(l1);

            ntripletsJGv=ntripletsJGv+1;
	    rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n2;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)-coefV_n2/(0.5*dt);

            rhs_G(l1) = rhs_G(l1)-sign_n2.*coef_n2'*A(ajlk_n2);
            rhs_G(l1) = rhs_G(l1)-coefV_n2*dtV(n2);        %if n1 n2 both inside,get dtV(n1)-dtV(n1),if only half,then dtV 
	    						   % and div A on surf satisfy lorentz gauge ,add up to 0
        end
							 %so the boundary effect of V come from this term and the following E1,dtE1

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
        for i = 1:length(ajvol_l1)
            switch ajvolM_l1(i)
                case 1
                    J   = (mJ_0(l1)+mJ_1(l1)+mJ_2(l1))+epsilon_mt*dtE1;
                    dJ_v = (dt*omega_p*omega_p/(dt*gamma_p+2)/linkL(l1)+...
			  (lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt)/linkL(l1) +...
			  (lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt)/linkL(l1)) +...
			  epsilon_mt/(0.5*dt*linkL(l1));
                    dJ_H = -(dt*omega_p*omega_p/(dt*gamma_p+2) +...
			   (lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt) +...
			   (lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt) +...
			   epsilon_mt/(0.5*dt));
                case 2
                    dJ_v = epsilon_in/(0.5*dt*linkL(l1));
                    dJ_H = -epsilon_in/(0.5*dt);
                    J = epsilon_in*dtE1+Js(l1);
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
     end
    
JG_v=sparse(rowJGv(1:ntripletsJGv),colJGv(1:ntripletsJGv),valJGv(1:ntripletsJGv),Nlink,Nnode);
JG_H=sparse(rowJGH(1:ntripletsJGH),colJGH(1:ntripletsJGH),valJGH(1:ntripletsJGH),Nlink,Nlink);
clear rowJGv colJGv valJGv rowJGH colJGH valJGH;
JG_v = JG_v(eqnLinks,eqnNodes);
JG_H  = JG_H(eqnLinks,eqnLinks);
JG = [JG_v,JG_H];
clear JG_v JG_H;
    
display(['time for matrix collection,G matrix:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rhs_F = rhs_F(eqnNodes);
    rhs_G = rhs_G(eqnLinks);
      
    Jacob =[JF;JG];
    rhs = -sparse([rhs_F;rhs_G]);

    display(['time for matrix reconstruction:']);
    toc;
    tic;

    clear JF  JG;
    clear rhs_F rhs_G;
%     Jacob1 = sparse(Jacob);
%     rhs1 = sparse(rhs);
%     [L,U] = luinc(Jacob1,droptol);
%     [dX,flag,relres,iter] = gmres(Jacob1,rhs1,[],linSolveTol,maxLinIt,L,U);
%     if flag ~= 0, error('No convergence for linear solver'); end
%     gmresItNr(l) = iter(2);
%    invJacob = Jacob^(-1);
%    dX = invJacob*rhs;
%    dX = Jacob\rhs; newton raphson: x_n+1 = x_n - f(x_n)/f'(x_n)
    dX = mylusovle(Jacob,rhs,1);
    display(['time for solving linear equation:']);
    toc;

    dV = dX(1:NeqnNodes);
    dH = dX(NeqnNodes+1:NeqnNodes+Nl);
    
    normRes = norm(rhs);
    normRes_pre = normRes;

    itNr = itNr+1;
    
    V(eqnNodes) = V(eqnNodes)+dV;
    H(eqnLinks) = H(eqnLinks)+dH;

    dtV  = 2 * (V - V_p)/dt - dtV1;
    dtH  = 2 * (H - H_p)/dt - dtH1;

    A(eqnLinks) = A_p(eqnLinks) + (H(eqnLinks) + H_p(eqnLinks))*dt/2;

     normUpdate = max([normRes]);%,norm(dV)/norm(V(eqnNodes)),norm(dH)/norm(H(eqnLinks))]);

     maxdV=max([norm(dV)/norm(V(eqnNodes))]);
     maxdH=max([norm(dH)/norm(H(eqnLinks))]);
     
     display([' maxdV: ',num2str(maxdV)]);
     display([' maxdH: ',num2str(maxdH)]);

    display(['  Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
        '; ','normUpdate:',num2str(normUpdate),'; ','normdV:',num2str(norm(dV))]);

    maxmJ0=max([mJ_0]);
    maxmJ1=max([mJ_1]);
    maxmJ2=max([mJ_2]);
    display(['max mJ0:',num2str(maxmJ0)]);
    display(['max mJ1:',num2str(maxmJ1)]);
    display(['max mJ2:',num2str(maxmJ2)]);
    clear dX dV dn dp dH;
   if (itNr - maxNewtonIt == 0) && normUpdate > updateTol
	display(['reach maxNewtonIt and still not converge']);
   end 
    
end    %( end of while) 
toc(tStart);   
display(['  End TD simulation relaxation2 per step.']);
