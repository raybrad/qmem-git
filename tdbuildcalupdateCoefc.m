function tdbuildcalupdateCoefc(dt)

global sigma epsilon_in;
global epsilon_mt;
global omega_p gamma_p;
global lepsr_1 lomega_1 lgamma_1 lepsr_2 lomega_2 lgamma_2;
global scl;
global links;
global Nnode Nlink;
global nodeLinks linkSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL volumeM;
global bndNodes dirNodes;
global isBndNodes;
global bndLinks;
global metalLinks;
global light_speed;
global EsurfLinks BsurfLinks;
%global Fn1matrix Fn2matrix Flkmatrix Fc;
global Fn1matrix;
%global Fc11;
global Fc;
global Gajlkmatrix11 Gajlk_n1matrix Gajlk_n2matrix Gn1matrix Gn2matrix;
%global Gc11 Gc12 Gc20 Gc31 Gc32 Gc41 Gc42 Gc51 Gc53 Gc54;
global Gc;

%%%%%%% initial guess  %%%%%%%%%%%%%


eqnNodes = setdiff((1:Nnode)',dirNodes);
eqnLinks=setdiff((1:Nlink)',EsurfLinks);	%for E boundary condition
isBsurfLinks=false(Nlink,1);
isBsurfLinks(BsurfLinks)=true;		        %for B boundary condition

%%%%%%%% start Newton's iteration %%%%%%%%
Eps = [epsilon_mt,epsilon_in,epsilon_qm];
prefac=[1,0,0];

NeqnNodes = length(eqnNodes);
Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
itNr = 0;

tStart=tic;
tic;
%%%% Lorezentz gauge%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% par               par     1                      % 
%---- V= div*A  --> ----V+ ---- integral0(A)dS=0   % 
%par t              par t   Vol                    %
%get (dV/dt)                                       %
%%%% Lorezentz gauge%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Fc11=zeros(Nnode,6);
Fn1matrix=zeros(Nnode,6);
for k = 1:NeqnNodes
    n1 = eqnNodes(k);
    ajlk_n1 = nodeLinks{n1}(1,:);
    ajnd_n1 = nodeLinks{n1}(2,:);
    ajvol_n1 = nodeVolumes{n1}(1,:);
    ajvolV_n1 = nodeVolumes{n1}(2,:);
    ajvolM_n1 = volumeM(ajvol_n1);
    sign_n1 = sign(ajnd_n1-n1);
    
    for i = 1:length(ajlk_n1)
      lk = ajlk_n1(i);
      coef_n1 = linkS(lk)/nodeV(n1);
      Fc11(n1,i)=-sign_n1(i)*coef_n1;
      Fn1matrix(n1,i)=lk;
      %rhs_F(n1)  = -sign_n1(i)*coef_n1*A(lk);
    end
    epsnode= sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1);
    Fc11(n1,:)=Fc11(n1,:)/(scl.K*epsnode);
    %dtV(n1) = rhs_F(n1) / (scl.K * epsnode);
end
Fc=[Fc11];
      %FAmatrix(eqnNodes,:)=A(Fn1matrix(eqnNodes,:));
      %rhs_F=sum((Fc11.*FAmatrix,2);
    
     % F= div J_tot = 0 (for node, no A) J_tot =J_f +epsilon part E/par t    
     % current continuity equation + gauss' law: div J_f+ d (div D)/dt = 0
     % => J_f* S+ eps *dE/dt *S =0
     % sigma*E + par E/par t (metal)                                 
     %for dielectric
     %

%for k = 1:NeqnNodes
%    n1 = eqnNodes(k);
%    ajlk_n1 = nodeLinks{n1}(1,:);
%    ajnd_n1 = nodeLinks{n1}(2,:);
%    ajvol_n1 = nodeVolumes{n1}(1,:);
%    ajvolV_n1 = nodeVolumes{n1}(2,:);
%    ajvolM_n1 = volumeM(ajvol_n1);
%    sign_n1 = sign(ajnd_n1-n1);
%    
%    for i = 1:length(ajlk_n1)
%        n2 = ajnd_n1(i);
%        lk = ajlk_n1(i);
%        ajvol_lk = linkVolumes{lk}(1,:);
%        ajvolS_lk = linkVolumes{lk}(2,:);
%        ajvolM_lk = volumeM(ajvol_lk);
%
%        dtE1 = -(dtVp(n2)-dtVp(n1))/linkL(lk)-sign_n1(i)*dtHp(lk);		
%                
%        epslink = sum(Eps(ajvolM_lk).*ajvolS_lk);
%
%       % rhs_F(n1)=rhs_F(n1)+sum(Eps(ajvolM_lk)*(dtVp(n2)/linkL(lk)+sign_n1(i)*dtHp(lk))*ajvolS_lk)/epslink*linkL(lk);
%		    
%        Fc21(n1,i)=sum(Eps(ajvolM_lk).*ajvolS_lk*(-1/linkL(lk))*(any(ajvolM_n1==1)));
%		Fc22(n1,i)=sum(Eps(ajvolM_lk).*ajvolS_lk*(-sign_n1(i))*(any(ajvolM_n1==1)));
%		
%        Fc30(n1,i)=sum(prefac(ajvolM_lk).*ajvolS_lk*sign_n1(i));			%mJ
%
%		Fn1matrix(n1,i)=n1;
%		Fn2matrix(n1,i)=n2;
%		Flkmatrix(n1,i)=lk;
%    end
%        %dtV(n1)=rhs_F(n1);
%        FdtVmatrix=-sum((Fc22.*FdtHmatrix+Fc30.*FJmatrix),2)./Fc21;
%end

Fc=[Fc11];
display(['time for matrix collection,F matrix:']);
toc;
tic;
%%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
Gc11=zeros(Nlink,12);
Gajlkmatrix11=ones(Nlink,12);

Gc12=zeros(Nlink,1);
Gc20=zeros(Nlink,1);

Gn1matrix=zeros(Nlink,1);
Gc31=zeros(Nlink,6);
Gajlk_n1matrix=ones(Nlink,6);
Gc32=zeros(Nlink,1);

Gn2matrix=zeros(Nlink,1);
Gc41=zeros(Nlink,6);
Gajlk_n2matrix=ones(Nlink,6);
Gc42=zeros(Nlink,1);

Gc51=zeros(Nlink,1);
%Gc52=zeros(Nlink,1);
Gc53=zeros(Nlink,1);
Gc54=zeros(Nlink,1);
Gn3matrix=zeros(Nlink,1);
Gc60=zeros(Nlink,6);
% -K*parH/part= cur curl A - K*{J_f+par E/par t} =0 (for link A,PI)   
 for k=1:Nl
    l1 = eqnLinks(k);
    n1 = links(l1,1);
    n2 = links(l1,2);
    n3 = links(l1,3);
    ajvol_l1 = linkVolumes{l1}(1,:);%volumeid
    ajvolS_l1 = linkVolumes{l1}(2,:);%associate surf area,not dual area, but only its own 1/4
    ajvolM_l1 = volumeM(ajvol_l1);
    CC = zeros(1,Nlink);

    ajlk_n1=[];
    ajlk_n2=[];
   
    Gn1matrix(l1)=n1;
    Gn2matrix(l1)=n2;
	Gn3matrix(l1)=n3;
    %%%%%%%%%% curl curl of A %%%%%%%%%%
    ajsf_l1 = linkSurfs{l1};
    Gc52_tmp=(isSQMlinks(l1)==0)*(-scl.K)*sum(Eps(ajvolM_l1).*ajvolS_l1*(-1))+...
	         (isSQMlinks(l1)~=0)*(-scl.K)*epsilon_qm*linkS(l1)*(-1);
    tmpcounter=0;
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
       % rhs_G(l1) = rhs_G(l1)+CC(ajlk)*A(ajlk);
    	Gajlkmatrix11(l1,tmpcounter+1:tmpcounter+3)=ajlk;
	    Gc11(l1,tmpcounter+1:tmpcounter+3)=CC(ajlk)/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);
        tmpcounter=tmpcounter+3;
    end
       %Gajlkmatrix11(l1,tmpcounter+1:12)=1;
       %Gc11(l1,:)=Gc11(l1,:).*(Gajlkmatrix11(l1,:)~=0);
       %rhs_G(l1) = rhs_G(l1)+CC(l1)*A(l1);
	   %rhs_G(l1) = rhs_G(l1)+((size(ajsf_l1,2)==3)*(1-isBsurfLinks(l1))*dL*(1/light_speed)+...
	   % 	     (size(ajsf_l1,2)==2)*(2-isBsurfLinks(l1))*dL*(1/light_speed))*H(l1);
	    Gc12(l1)=CC(l1)/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);
		
	    Gc20(l1)=((size(ajsf_l1,2)==3)*(1-isBsurfLinks(l1))*dL*(1/light_speed)+...
	    	     (size(ajsf_l1,2)==2)*(2-isBsurfLinks(l1))*dL*(1/light_speed))/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);
      if ~isBndNodes(n1)
          ajlk_n1 = nodeLinks{n1}(1,:); % sta node connected to at most 6 links
          ajnd_n1 = nodeLinks{n1}(2,:); % sta node connected to at most 6 nodes
          sign_n1 = sign(ajnd_n1-n1);
          coef_n1 = linkS(ajlk_n1)./nodeV(n1)*linkS(l1)/linkL(l1); %linkS(Nlink,1)~ dual surface area,nodeV(Nnode,1)
          ajvol_n1 = nodeVolumes{n1}(1,:);
          ajvolV_n1 = nodeVolumes{n1}(2,:);
          ajvolM_n1 = volumeM(ajvol_n1);
          coefV_n1 = scl.K*(sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1))*linkS(l1)/linkL(l1);
	
          %rhs_G(l1) = rhs_G(l1)-(-sign_n1.*coef_n1')*A(ajlk_n1); %gradient divergence A,if it's connected to a bndNode,linkS
	  	  %						   %could be on bndsurf
          %rhs_G(l1) = rhs_G(l1)+coefV_n1*dtVp(n1);		   %dV/dt
    	  Gc31(l1,1:length(ajlk_n1))=-(-sign_n1.*coef_n1')/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);
	      Gajlk_n1matrix(l1,1:length(ajlk_n1))=ajlk_n1;
	      Gc32(l1)=coefV_n1/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);
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

          %rhs_G(l1) = rhs_G(l1)-sign_n2.*coef_n2'*A(ajlk_n2);
          %rhs_G(l1) = rhs_G(l1)-coefV_n2*dtVp(n2);        %if n1 n2 both inside,get dtV(n1)-dtV(n1),if only half,then dtV 
	  	  %					   % and div A on surf satisfy lorentz gauge ,add up to 0
    	  Gc41(l1,1:length(ajlk_n2))=-(sign_n2.*coef_n2')/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);
    	  Gajlk_n2matrix(l1,1:length(ajlk_n2))=ajlk_n2;
    	  Gc42(l1)= -coefV_n2/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);
      end

    %%%%%% source current %%%%%%%%%%%
          epslink = sum(Eps(ajvolM_l1).*ajvolS_l1);
          %dtH(l1) =-(rhs_G(l1)-scl.K*sum(Eps(ajvolM_l1).*ajvolS_l1*(-1/linkL(l1)))*(dtVp(n2)-dtVp(n1))-scl.K*sum(ajvolS_l1)*Js(l1))/(-scl.K*sum(Eps(ajvolM_l1).*ajvolS_l1*(-1)));
       	  
       	  Gc51(l1)=-scl.K*((isSQMlinks(l1)==0)*sum(Eps(ajvolM_l1).*ajvolS_l1*(-1/linkL(l1)))+...
                           (isSQMlinks(l1)~=0)*epsilon_qm*linkS(l1)*(-1/linkL(l1)))/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);
       	  %Gc52(l1)= Gc52_tmp+(Gc52_tmp==0);
    	  Gc53(l1)=-scl.K*(isSQMlinks(l1)==0)*sum(prefac(ajvolM_l1).*ajvolS_l1)/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);	%mJ
    	  Gc54(l1)=-scl.K*(isSQMlinks(l1)==0)*sum(ajvolS_l1)/(Gc52_tmp+(Gc52_tmp==0))*(Gc52_tmp~=0);				%Js
	      Gc60(l1,n3+3)=(isSQMlinks(l1)~=0)*(-scl.K)*linkS(l1);				%currdlink
          %Gc Gc51*(dtVp(n2)-dtVp(n1))+Gc54*Js)/Gc52;
          %dtH =-sum(Gc11.*GAmatrix11+Gc12.*A+Gc20.*H+sum(Gc31.*GAmatrix31,2)+Gc32.*GdtVmatrix32+sum(Gc41.*GAmatrix41,2)+Gc42.*GdtVmatrix42+Gc51.*GdtVmatrix51+Gc53.*(mJ_0+mJ_1+mJ_2)+Gc54.*Js,2)/Gc52;  
 end

Gc=[Gc11,Gc12,Gc20,Gc31,Gc32,Gc41,Gc42,Gc51,Gc53,Gc54,Gc60];
display(['time for matrix collection,G matrix:']);
toc;

toc(tStart);   

display(['  End build tdcalupdatec Coefficient.']);
