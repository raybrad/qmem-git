function tdbuildJacobc(dt)

display(['  Start building Jacobc matrix:']);

global sigma epsilon_in;
global epsilon_mt;
global epsilon_qm;
global omega_p gamma_p;
global lepsr_1 lomega_1 lgamma_1 lepsr_2 lomega_2 lgamma_2;
global scl;
global links;
global Nnode Nlink;
global nodeLinks linkSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL volumeM;
global dirNodes;
global isBndNodes;
global metalLinks;

global isSQMlinks;



global light_speed;
global EsurfLinks BsurfLinks;
%global Jacob colind Lmatrix Umatrix;
global JacobLU;
%Method 2
%Solving Gauss' law.
%Solving Current-continuity equations.
%Solving Maxwell-Ampere equations(already insert Lorentz gauge)
%Respecting the gauge condition is a side product of solving above set.
%%%%%%% initial guess  %%%%%%%%%%%%%
eqnNodes = setdiff((1:Nnode)',dirNodes);%eqnNodes (nodes except from contact ,different from tdcalupdatec.m)   
					%dirNodes (contact part update by dtVp)

eqnLinks=setdiff((1:Nlink)',EsurfLinks);	%for E boundary condition
isBsurfLinks=false(Nlink,1);
isBsurfLinks(BsurfLinks)=true;		        %for B boundary condition

%%%%%%%% start Newton's iteration %%%%%%%%
updateTol = 1e-4;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 50;
maxLinIt = 200;
Eps = [epsilon_mt,epsilon_in,epsilon_qm];

NeqnNodes = length(eqnNodes);
Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
itNr = 0;
tStart=tic;

%JF
% construct the sparse2 matrix by (Row,Column,Value) 
tic;
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for k=1:NeqnNodes
        n1 = eqnNodes(k);
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
        sign_n1 = sign(ajnd_n1-n1);

        if all(ajvolM_n1 == 3)	%all nodes in QM region (include the links along the boundary)
            for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    n3 = links(lk,3);
               %    dtE1 = -(dtV(n2)-dtV(n1))/linkL(lk)-sign_n1(i)*dtH(lk);
               %    rhs_F(n1)   = rhs_F(n1) + sign_n1(i)*linkS(lk)*currdlink(lk,n3+3) + linkS(lk)* epsilon_qm *dtE1;
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
        elseif any(ajvolM_n1 == 1) % metal/semi interfaces or triple points
            for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    n3 = links(lk,3);
%                   dtE1 = -(dtV(n2)-dtV(n1))/linkL(lk)-sign_n1(i)*dtH(lk);
                    if isSQMlinks(lk) == 0	
	            %along(impossible here) or connected to the outside QM boundary,or pure EM links 
                      ajvol_lk = linkVolumes{lk}(1,:);
                      ajvolS_lk = linkVolumes{lk}(2,:);
                      ajvolM_lk = volumeM(ajvol_lk);
 
                      for j = 1:length(ajvol_lk)
                         switch ajvolM_lk(j)
                            case 1
                  		%J=sign_n1(i)*(mJ_0(lk)+mJ_1(lk)+mJ_2(lk))+epsilon_mt*dtE1;
				dJv =(dt*omega_p*omega_p/(dt*gamma_p+2)/linkL(lk)+...
				    (lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt)/linkL(lk) +...
				    (lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt)/linkL(lk)) +...
				     epsilon_mt/(0.5*dt*linkL(lk));
                                dJH = -sign_n1(i)*(dt*omega_p*omega_p/(dt*gamma_p+2) +...
				                  (lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt) +...
				                  (lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt) +...
						  epsilon_mt/(0.5*dt));
                            case 2
%                                J   = epsilon_in*dtE1;
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

%                        rhs_F(n1) = rhs_F(n1)+ajvolS_lk(j)*J;
                      end
                    elseif isSQMlinks(lk) <= 6	%inside/outside interface
%                      rhs_F(n1)   = rhs_F(n1) + sign_n1(i)*linkS(lk)*currdlink(lk,n3+3) + linkS(lk)* epsilon_qm *dtE1;
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
%                    elseif isSQMlinks(lk) > 6 %inside QM, very rare in ajvolM_n1 == 2 not ajvolM_n1==4 ,it's 0 anyway
%                       rhs_F(n1)   = rhs_F(n1) + 0;
                    end
            end
        else   %other ajvolM_n1 except QMnodes and metal 
                for i = 1:length(ajlk_n1)
                    n2 = ajnd_n1(i);
                    lk = ajlk_n1(i);
                    if isSQMlinks(lk) == 0  %along(impossible here)or connected to the outside QM boundary,or pure EM links 
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
                      coefH = -coefV*linkL(lk);
		      ntripletsJFH=ntripletsJFH+1;
		      rowJFH(ntripletsJFH)=n1;
		      colJFH(ntripletsJFH)=lk;
		      valJFH(ntripletsJFH)=sign_n1(i)*coefH;
%rhs_F(n1) = rhs_F(n1)+sum(Eps(ajvolM_lk).*(-(V(n2)-V(n1))/linkL(lk)-sign_n1(i)* H(lk)).*ajvolS_lk);
                    else	%interface
%                     rhs_F(n1)   = rhs_F(n1) + epsilon_qm *(-(V(n2)-V(n1))/linkL(lk)-sign_n1(i)* H(lk))*linkS(lk);
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
        end
     end
JF_v=sparse2(rowJFv(1:ntripletsJFv),colJFv(1:ntripletsJFv),valJFv(1:ntripletsJFv),Nnode,Nnode);
JF_H=sparse2(rowJFH(1:ntripletsJFH),colJFH(1:ntripletsJFH),valJFH(1:ntripletsJFH),Nnode,Nlink);
clear rowJFv colJFv valJFv rowJFH colJFH valJFH;
JF_v = JF_v(eqnNodes,eqnNodes);
JF_H = JF_H(eqnNodes,eqnLinks);
JF = [JF_v,JF_H];
clear JF_v JF_H;
display(['time for matrix collection,JF matrix:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        n1 = links(l1,1);
        n2 = links(l1,2);
        ajvol_l1 = linkVolumes{l1}(1,:);%volumeid
        ajvolS_l1 = linkVolumes{l1}(2,:);%associate surf area,not dual area, but only its own 1/4
        ajvolM_l1 = volumeM(ajvol_l1);
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
            %rhs_G(l1) = rhs_G(l1)+CC(ajlk)*A(ajlk);
	    for m=1:length(ajlk)
	    ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=ajlk(m);
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+CC(ajlk(m))*dt/2;
	    end
        end

        %rhs_G(l1) = rhs_G(l1)+CC(l1)*A(l1);
	ntripletsJGH=ntripletsJGH+1;
	rowJGH(ntripletsJGH)=l1;
	colJGH(ntripletsJGH)=l1;
	valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+CC(l1)*dt/2;

	    if (size(ajsf_l1,2)==3) %xy plane : Ay(lk)~Bx~par Ay(lk)/par z =1/c par Ay(lk)/par t no matter top or bottom
	        if(isBsurfLinks(l1))
%		rhs_G(l1)=rhs_G(l1)+0.0;
		bndJG_H(l1)=0.0;
	        ntripletsJGH=ntripletsJGH+1;
		rowJGH(ntripletsJGH)=l1;
		colJGH(ntripletsJGH)=l1;
		valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
       	    else
%       		rhs_G(l1)=rhs_G(l1)+dL*(1/light_speed*H(l1));
       		bndJG_H(l1)=dL/light_speed;
	        ntripletsJGH=ntripletsJGH+1;
		rowJGH(ntripletsJGH)=l1;
		colJGH(ntripletsJGH)=l1;
		valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
	    end
	    elseif (size(ajsf_l1,2)==2) 
	        if(isBsurfLinks(l1))
%		rhs_G(l1)=rhs_G(l1)+dL*(1/light_speed*H(l1));
		bndJG_H(l1)=dL/light_speed;
	        ntripletsJGH=ntripletsJGH+1;
		rowJGH(ntripletsJGH)=l1;
		colJGH(ntripletsJGH)=l1;
		valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
		else
%		rhs_G(l1)=rhs_G(l1)+2*dL*(1/light_speed*H(l1));
		bndJG_H(l1)=2*dL/light_speed;
	        ntripletsJGH=ntripletsJGH+1;
		rowJGH(ntripletsJGH)=l1;
		colJGH(ntripletsJGH)=l1;
		valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+bndJG_H(l1);
		end
	    end
        %%%%%%%% gradient divergence of A and gradient of V %%%%%%%%
        if ~isBndNodes(n1) 
            ajlk_n1 = nodeLinks{n1}(1,:); % sta node connected to at most 6 links
            ajnd_n1 = nodeLinks{n1}(2,:); % sta node connected to at most 6 nodes
            sign_n1 = sign(ajnd_n1-n1);
            coef_n1 = linkS(ajlk_n1)./nodeV(n1)*linkS(l1)/linkL(l1);
            GD(ajlk_n1) = GD(ajlk_n1)-sign_n1.*coef_n1';
            ajvol_n1 = nodeVolumes{n1}(1,:);
            ajvolV_n1 = nodeVolumes{n1}(2,:);
            ajvolM_n1 = volumeM(ajvol_n1);
            coefV_n1 = scl.K*(sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1))*linkS(l1)/linkL(l1);

            ntripletsJGv=ntripletsJGv+1;
	    rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n1;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+coefV_n1/(0.5*dt);

            %rhs_G(l1) = rhs_G(l1)-(-sign_n1.*coef_n1')*A(ajlk_n1);
            %rhs_G(l1) = rhs_G(l1)+coefV_n1*dtV(n1);
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

            %rhs_G(l1) = rhs_G(l1)-sign_n2.*coef_n2'*A(ajlk_n2);
            %rhs_G(l1) = rhs_G(l1)-coefV_n2*dtV(n2);
        end

	ajlk_n1n2=union(ajlk_n1,ajlk_n2);
        for m=1:length(ajlk_n1n2)
	ntripletsJGH=ntripletsJGH+1;
	rowJGH(ntripletsJGH)=l1;
	colJGH(ntripletsJGH)=ajlk_n1n2(m);
	valJGH(ntripletsJGH)=valJGH(ntripletsJGH)-GD(ajlk_n1n2(m))*dt/2;
	end

        %%%%%% source current %%%%%%%%%%%
       % E1 = -(V(n2)-V(n1))/linkL(l1)-H(l1); %(div V+ PI)
       % dtE1 = -(dtV(n2)-dtV(n1))/linkL(l1) - dtH(l1);
        if isSQMlinks(l1) == 0		 %along or connected to the outside QM boundary,or pure EM links
       % E1 = -(V(n2)-V(n1))/linkL(l1)-H(l1);
          for i = 1:length(ajvol_l1)	%volumetype of the link belong to
            switch ajvolM_l1(i)
                case 1
                    %J   = (mJ_0(l1)+mJ_1(l1)+mJ_2(l1))+epsilon_mt*dtE1;
                    dJv = (dt*omega_p*omega_p/(dt*gamma_p+2)/linkL(l1)+...
			  (lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt)/linkL(l1) +...
			  (lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt)/linkL(l1)) +...
			  epsilon_mt/(0.5*dt*linkL(l1));
                    dJH = -(dt*omega_p*omega_p/(dt*gamma_p+2) +...
			   (lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt) +...
			   (lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt) +...
			   epsilon_mt/(0.5*dt));
                case 2
                    %J = epsilon_in*dtE1+Js(l1);
                    dJv = epsilon_in/(0.5*dt*linkL(l1));
                    dJH = -epsilon_in/(0.5*dt);
                case 3	%all QM nodes(include the outside/inside boundary)
                    %J = epsilon_qm*dtE1;
                    dJv = epsilon_qm/(0.5*dt*linkL(l1));
                    dJH = -dJv*linkL(l1);
                otherwise
                    error('undefined material');
            end
            ntripletsJGv=ntripletsJGv+1;
            rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n1;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)-scl.K*ajvolS_l1(i)*dJv;
	    ntripletsJGv=ntripletsJGv+1;
	    rowJGv(ntripletsJGv)=l1;
	    colJGv(ntripletsJGv)=n2;
	    valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+scl.K*ajvolS_l1(i)*dJv;

            ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=l1;
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)-scl.K*ajvolS_l1(i)*dJH;
%            rhs_G(l1) = rhs_G(l1)-scl.K*ajvolS_l1(i)*J;
	  end
        else %mix  outside/inside QM links interface, inside QM links
          %  rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*(currdlink(l1,n3+3)+epsilon_qm*dtE1);
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
    
JG_v=sparse2(rowJGv(1:ntripletsJGv),colJGv(1:ntripletsJGv),valJGv(1:ntripletsJGv),Nlink,Nnode);
JG_H=sparse2(rowJGH(1:ntripletsJGH),colJGH(1:ntripletsJGH),valJGH(1:ntripletsJGH),Nlink,Nlink);
clear rowJGv colJGv valJGv rowJGH colJGH valJGH;
JG_v = JG_v(eqnLinks,eqnNodes);
JG_H  = JG_H(eqnLinks,eqnLinks);
JG = [JG_v,JG_H];
clear JG_v JG_H;
    
display(['time for matrix collection,JG matrix:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jacob =sparse2([JF;JG]);

clear JF JG;

%colind=colamd(Jacob);
%aspas = sparse2(Jacob(:,colind));
%[Lmatrix,Umatrix] = ilu(aspas,struct('type','ilutp','droptol',1e-4));	%most time consuming
%savefilename='JacobMatrix.mat';
%save(savefilename, 'Jacob','colind','Lmatrix', 'Umatrix','-v7.3');

[JacobLU.L, JacobLU.U, JacobLU.P, JacobLU.Q, JacobLU.R] = lu (Jacob) ;
%savefilename='JacobLU.mat';
%save(savefilename, 'JacobLU','-v7.3');

display(['time for matrix Jacob matrix combination:']);
toc;

toc(tStart);   
display(['  End building Jacob matrix.']);
