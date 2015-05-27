function tdbuildRHSCoef(dt)

display(['  Start building RHS coefficientsc matrix:']);

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

%global Fc11 Fc12 Fc21 Fc22 Fc30 Fn1matrix Fn2matrix Flkmatrix;
%global Gc11 Gajlkmatrix11 Gc12 Gc20 Gn1matrix Gc31 Gajlk_n1matrix Gc32;
%global Gn2matrix Gc41 Gajlk_n2matrix Gc42 Gc51 Gc52 Gc53 Gc54;
global Fn1matrix Fn2matrix Fn3matrix Flkmatrix Fc;
global Gajlkmatrix11 Gajlk_n1matrix Gajlk_n2matrix Gn1matrix Gn2matrix Gn3matrix Gc;
%Method 2
%Solving Gauss' law.
%Solving Current-continuity equations.
%Solving Maxwell-Ampere equations(already insert Lorentz gauge)
%Respecting the gauge condition is a side product of solving above set.

eqnNodes = setdiff((1:Nnode)',dirNodes);%eqnNodes (nodes except from contact ,different from tdcalupdatec.m)   
					%dirNodes (contact part update by dtVp)

eqnLinks=setdiff((1:Nlink)',EsurfLinks);	%for E boundary condition
isBsurfLinks=false(Nlink,1);
isBsurfLinks(BsurfLinks)=true;		        %for B boundary condition

%%%%%%%% start Newton's iteration %%%%%%%%
Eps = [epsilon_mt,epsilon_in,epsilon_qm];
prefac=[1,0,0];

NeqnNodes = length(eqnNodes);
Nl = length(eqnLinks);

tStart=tic;

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rhsF
% construct the sparse2 matrix by (Row,Column,Value) 
 rhs_F = zeros(Nnode,1);
 Fc11=zeros(Nnode,6);
 Fc12=zeros(Nnode,6);
 Fc21=zeros(Nnode,6);
 Fc22=zeros(Nnode,6);
 Fc30=zeros(Nnode,6);
 Fc40=zeros(Nnode,6);
 Fn1matrix =zeros(Nnode,6);
 Fn2matrix =zeros(Nnode,6);
 Fn3matrix =zeros(Nnode,6);
 Flkmatrix =zeros(Nnode,6);

for k=1:NeqnNodes
        n1 = eqnNodes(k);
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
                    ajvol_lk = linkVolumes{lk}(1,:);
                    ajvolS_lk = linkVolumes{lk}(2,:);
                    ajvolM_lk = volumeM(ajvol_lk);

		    Fc11(n1,i)=(all(ajvolM_n1==3))*epsilon_qm*linkS(lk)*(-1/linkL(lk))+ ...
		               (any(ajvolM_n1==1))*(sum((isSQMlinks(lk)==0)*Eps(ajvolM_lk).*ajvolS_lk*(-1/linkL(lk)))+...
		                                        (isSQMlinks(lk)>=1 && isSQMlinks(lk)<=6)*epsilon_qm*linkS(lk)*(-1/linkL(lk)));

		    Fc12(n1,i)=(all(ajvolM_n1==3))*epsilon_qm*linkS(lk)*(-sign_n1(i))+...
		               (any(ajvolM_n1==1))*(sum((isSQMlinks(lk)==0)*Eps(ajvolM_lk).*ajvolS_lk*(-sign_n1(i))) +...
		                                        (isSQMlinks(lk)>=1 && isSQMlinks(lk)<=6)*epsilon_qm*linkS(lk)*(-sign_n1(i)));

		    Fc21(n1,i)=(any(ajvolM_n1~=3) && all(ajvolM_n1~=1))*(sum((isSQMlinks(lk)==0)*Eps(ajvolM_lk).*ajvolS_lk*(-1/linkL(lk)))+...
		                   					     (isSQMlinks(lk)~=0)*epsilon_qm*linkS(lk)*(-1/linkL(lk)));

		    Fc22(n1,i)=(any(ajvolM_n1~=3) && all(ajvolM_n1~=1))*(sum((isSQMlinks(lk)==0)*Eps(ajvolM_lk).*ajvolS_lk*(-sign_n1(i))) +... 	
		                    				             (isSQMlinks(lk)~=0)*epsilon_qm*linkS(lk)*(-sign_n1(i)));

		    Fc30(n1,i)=(any(ajvolM_n1==1) && isSQMlinks(lk)==0)*sum(prefac(ajvolM_lk).*ajvolS_lk*sign_n1(i));	%mJ

	            Fc40(n1,n3+3)=(all(ajvolM_n1==3)+((any(ajvolM_n1==1) && isSQMlinks(lk)>=1 && isSQMlinks(lk)<=6)))*sign_n1(i)*linkS(lk); %currdlink

		    Fn1matrix(n1,i)=n1;
		    Fn2matrix(n1,i)=n2;
		    Fn3matrix(n1,i)=n3;
		    Flkmatrix(n1,i)=lk;
            end
end
display(['time for rhs_F coefficient matrix collection:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rhsG
rhs_G = zeros(Nlink,1);
Gc11=zeros(Nlink,12);
Gajlkmatrix11=zeros(Nlink,12);

Gc12=zeros(Nlink,1);
Gc20=zeros(Nlink,1);

Gn1matrix=zeros(Nlink,1);
Gc31=zeros(Nlink,6);
Gajlk_n1matrix=zeros(Nlink,6);
Gc32=zeros(Nlink,1);

Gn2matrix=zeros(Nlink,1);
Gc41=zeros(Nlink,6);
Gajlk_n2matrix=zeros(Nlink,6);
Gc42=zeros(Nlink,1);

Gc51=zeros(Nlink,1);
Gc52=zeros(Nlink,1);
Gc53=zeros(Nlink,1);
Gc54=zeros(Nlink,1);
Gn3matrix=zeros(Nlink,1);
Gc60=zeros(Nlink,6);


%%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
    % G = cur curl A - {grad div A+grad V}+{-J+par E/par t} =0 (for link A,PI)
    %  already contains lorentz gauge
    %PI_n+1 =PI_n - G/{par G/parPI}  to get the PI

    for k=1:Nl 
        l1 = eqnLinks(k);
        n1 = links(l1,1);
        n2 = links(l1,2);
        ajvol_l1 = linkVolumes{l1}(1,:);%volumeid
        ajvolS_l1 = linkVolumes{l1}(2,:);%associate surf area,not dual area, but only its own 1/4
        ajvolM_l1 = volumeM(ajvol_l1);
        n3 = links(l1,3);
        CC = zeros(1,Nlink);
	ajlk_n1=[];
	ajlk_n2=[];
	Gn1matrix(l1)=n1;
	Gn2matrix(l1)=n2;
	Gn3matrix(l1)=n3;
        %%%%%%%%%% curl curl of A %%%%%%%%%%
        ajsf_l1 = linkSurfs{l1};
	tmpcounter=0;
        for i = 1:size(ajsf_l1,2)
            ajlk = ajsf_l1(2:4,i);
            S = linkL(l1)*linkL(ajlk(1));
            dL = dlinkL(ajsf_l1(1,i));
            ajlkL = linkL(ajlk)';
            
	    ajlkSign = (n2==links(ajlk(1),1))*[1,-1,-1]+(n2==links(ajlk(1),2))*[-1,-1,1];

            CC(ajlk) = dL/S*(ajlkSign.*ajlkL);
            CC(l1) = CC(l1)+dL/S*linkL(l1);
	    
	    Gajlkmatrix11(l1,tmpcounter+1:tmpcounter+3)=ajlk;
	    Gc11(l1,tmpcounter+1:tmpcounter+3)=CC(ajlk);
	    tmpcounter=tmpcounter+3;
            %rhs_G(l1) = rhs_G(l1)+CC(ajlk)*A(ajlk);
        end
            %rhs_G(l1) = rhs_G(l1)+CC(l1)*A(l1);
	    Gc12(l1)=CC(l1);
		
	    Gc20(l1)=(size(ajsf_l1,2)==3)*(1-isBsurfLinks(l1))*dL*(1/light_speed)+...
	    	     (size(ajsf_l1,2)==2)*(2-isBsurfLinks(l1))*dL*(1/light_speed);


        %%%%%%%% gradient divergence of A and gradient of V %%%%%%%%
        if ~isBndNodes(n1) 
            ajlk_n1 = nodeLinks{n1}(1,:); % sta node connected to at most 6 links
            ajnd_n1 = nodeLinks{n1}(2,:); % sta node connected to at most 6 nodes
            sign_n1 = sign(ajnd_n1-n1);
            coef_n1 = linkS(ajlk_n1)./nodeV(n1)*linkS(l1)/linkL(l1);
            ajvol_n1 = nodeVolumes{n1}(1,:);
            ajvolV_n1 = nodeVolumes{n1}(2,:);
            ajvolM_n1 = volumeM(ajvol_n1);
            coefV_n1 = scl.K*(sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1))*linkS(l1)/linkL(l1);
	
	    Gc31(l1,1:length(ajlk_n1))=-(-sign_n1.*coef_n1');
	    Gajlk_n1matrix(l1,1:length(ajlk_n1))=ajlk_n1;
	    Gc32(l1)=coefV_n1;
            %rhs_G(l1) = rhs_G(l1)-(-sign_n1.*coef_n1')*A(ajlk_n1);
            %rhs_G(l1) = rhs_G(l1)+coefV_n1*dtV(n1);
        end

        if ~isBndNodes(n2) 
            ajlk_n2 = nodeLinks{n2}(1,:);%6 neighbor links around node n2
            ajnd_n2 = nodeLinks{n2}(2,:);%6 neighboring nodes around node n2 
            sign_n2 = sign(ajnd_n2-n2);
            coef_n2 = linkS(ajlk_n2)./nodeV(n2)*linkS(l1)/linkL(l1);
            ajvol_n2 = nodeVolumes{n2}(1,:);	%volumeid
            ajvolV_n2 = nodeVolumes{n2}(2,:);	%associate 1/8 volume
            ajvolM_n2 = volumeM(ajvol_n2);
            coefV_n2 = scl.K*(sum(Eps(ajvolM_n2).*ajvolV_n2)/nodeV(n2))*linkS(l1)/linkL(l1);

	    Gc41(l1,1:length(ajlk_n2))=-(sign_n2.*coef_n2');
	    Gajlk_n2matrix(l1,1:length(ajlk_n2))=ajlk_n2;
	    Gc42(l1)= -coefV_n2;
            %rhs_G(l1) = rhs_G(l1)-sign_n2.*coef_n2'*A(ajlk_n2);
            %rhs_G(l1) = rhs_G(l1)-coefV_n2*dtV(n2);
        end

       	Gc51(l1)=(isSQMlinks(l1)==0)*(-scl.K)*sum(Eps(ajvolM_l1).*ajvolS_l1*(-1/linkL(l1)))+...
	         (isSQMlinks(l1)~=0)*(-scl.K)*epsilon_qm*linkS(l1)*(-1/linkL(l1));
       	Gc52(l1)=(isSQMlinks(l1)==0)*(-scl.K)*sum(Eps(ajvolM_l1).*ajvolS_l1*(-1))+...
	         (isSQMlinks(l1)~=0)*(-scl.K)*epsilon_qm*linkS(l1)*(-1);
	Gc53(l1)=(isSQMlinks(l1)==0)*(-scl.K)*sum(prefac(ajvolM_l1).*ajvolS_l1);		%mJ
	Gc54(l1)=(isSQMlinks(l1)==0)*(-scl.K)*sum(ajvolS_l1);							%Js
	Gc60(l1,n3+3)=(isSQMlinks(l1)~=0)*(-scl.K)*linkS(l1);				%currdlink
     end
    
display(['time for rhs_G coefficient matrix collection:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fc=[Fc11,Fc12,Fc21,Fc22,Fc30,Fc40];

Gc=[Gc11,Gc12,Gc20,Gc31,Gc32,Gc41,Gc42,Gc51,Gc52,Gc53,Gc54,Gc60];

display(['time for saving matrix: ']);
toc;
tic;

toc(tStart);   
display(['  End RHS coefficient matrix collection.']);
