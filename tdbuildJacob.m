function tdbuildJacob(dt)

display(['  Start building Jacob matrix:']);

global epsilon_in;
global epsilon_mt;
%
global scl;
global CJ01 CJ02 CJ11 CJ12 CJ13 CJ21 CJ22 CJ23;
global links;
global Nnode Nlink;
global nodeLinks linkSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL volumeM;
global bndNodes dirNodes;
global isBndNodes;
global bndLinks;
global metalLinks;
global light_speed XZsurfLinks XYsurfLinks YZsurfLinks;
global EsurfLinks BsurfLinks;
%global Jacob colind Lmatrix Umatrix;
global JacobLU permutation;
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
updateTol = 1e-8;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 50;
maxLinIt = 200;
Eps = [epsilon_mt,epsilon_in];
prefac=[1,0];

NeqnNodes = length(eqnNodes);
Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
itNr = 0;
tStart=tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                    ajvol_lk = linkVolumes{lk}(1,:);
                    ajvolS_lk = linkVolumes{lk}(2,:);
                    ajvolM_lk = volumeM(ajvol_lk);
				
			
			ntripletsJFv=ntripletsJFv+1;
			rowJFv(ntripletsJFv)=n1;
			colJFv(ntripletsJFv)=n1;
			valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+ ...
			     		     any(ajvolM_n1==1)*sum((prefac(ajvolM_lk)*(CJ02+CJ13+CJ23)+Eps(ajvolM_lk)/(0.5*dt)).*ajvolS_lk)/linkL(lk)+...
				             all(ajvolM_n1~=1)*sum(Eps(ajvolM_lk).*ajvolS_lk)/linkL(lk);
			ntripletsJFv=ntripletsJFv+1;
			rowJFv(ntripletsJFv)=n1;
			colJFv(ntripletsJFv)=n2;
			valJFv(ntripletsJFv)=valJFv(ntripletsJFv)- ...
			     		     any(ajvolM_n1==1)*sum((prefac(ajvolM_lk)*(CJ02+CJ13+CJ23)+Eps(ajvolM_lk)/(0.5*dt)).*ajvolS_lk)/linkL(lk)-...
				             all(ajvolM_n1~=1)*sum(Eps(ajvolM_lk).*ajvolS_lk)/linkL(lk);

			ntripletsJFH=ntripletsJFH+1;
			rowJFH(ntripletsJFH)=n1;
			colJFH(ntripletsJFH)=lk;
			valJFH(ntripletsJFH)=valJFH(ntripletsJFH)+ ...
			     		     any(ajvolM_n1==1)*(-sign_n1(i))*sum((prefac(ajvolM_lk)*(CJ02+CJ13+CJ23)+Eps(ajvolM_lk)/(0.5*dt)).*ajvolS_lk)+...
				             all(ajvolM_n1~=1)*(-sign_n1(i))*sum(Eps(ajvolM_lk).*ajvolS_lk);

            end
     end
JF_v=sparse2(rowJFv(1:ntripletsJFv),colJFv(1:ntripletsJFv),valJFv(1:ntripletsJFv),Nnode,Nnode);
JF_H=sparse2(rowJFH(1:ntripletsJFH),colJFH(1:ntripletsJFH),valJFH(1:ntripletsJFH),Nnode,Nlink);
clear rowJFv colJFv valJFv rowJFH colJFH valJFH;
JF_v = JF_v(eqnNodes,eqnNodes);
JF_H = JF_H(eqnNodes,eqnLinks);
%savefilename='JF.mat';
%save(savefilename, 'JF_v','JF_H','-v7.3');
JF = [JF_v,JF_H];
clear JF_v JF_H;
display(['time for matrix collection,JF matrix:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%JG
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
        for i = 1:size(ajsf_l1,2)	  %num of cell ,i.e. num of surfs the link is on
            ajlk = ajsf_l1(2:4,i);	  % the other 3 link on that surface i
            S = linkL(l1)*linkL(ajlk(1)); %surface area the link l1 is on, linkL ~link length  
            dL = dlinkL(ajsf_l1(1,i));    %1~surf index,dLinkL~dual length of links, (Nsurf,1)
            ajlkL = linkL(ajlk)';	  % length of the other 3 link
            if n2 == links(ajlk(1),1)	  %the end node of link l1 is the start node of the first surround link 
                ajlkSign = [1,-1,-1];
            elseif n2 == links(ajlk(1),2)
                ajlkSign = [-1,-1,1];
            else
                error('Incorrect link arrangement');
            end


            CC(ajlk) = dL/S*(ajlkSign.*ajlkL); %CC(1:3) other 3 link ,those link could be on the boundary,just l1 is not bnd
            CC(l1) = CC(l1)+dL/S*linkL(l1);	% link l1
	    for m=1:length(ajlk)
	    ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=ajlk(m);
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+CC(ajlk(m))*dt/2;
	    end
        end

	ntripletsJGH=ntripletsJGH+1;
	rowJGH(ntripletsJGH)=l1;
	colJGH(ntripletsJGH)=l1;
	valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+CC(l1)*dt/2;
%
        ntripletsJGH=ntripletsJGH+1;
        rowJGH(ntripletsJGH)=l1;
        colJGH(ntripletsJGH)=l1;
        valJGH(ntripletsJGH)=valJGH(ntripletsJGH)+ ...
      			    (size(ajsf_l1,2)==3)*(1-isBsurfLinks(l1))*dL*(1/light_speed)+...
         	    	    (size(ajsf_l1,2)==2)*(2-isBsurfLinks(l1))*dL*(1/light_speed);



        
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
        end
							 %so the boundary effect of V come from this term and the following E1,dtE1

	    ajlk_n1n2=union(ajlk_n1,ajlk_n2);
            for m=1:length(ajlk_n1n2)
	    ntripletsJGH=ntripletsJGH+1;
	    rowJGH(ntripletsJGH)=l1;
	    colJGH(ntripletsJGH)=ajlk_n1n2(m);
	    valJGH(ntripletsJGH)=valJGH(ntripletsJGH)-GD(ajlk_n1n2(m))*dt/2;
	    end
	
        %%%%%% source current %%%%%%%%%%%
	ntripletsJGv=ntripletsJGv+1;
	rowJGv(ntripletsJGv)=l1;
	colJGv(ntripletsJGv)=n1;
	valJGv(ntripletsJGv)=valJGv(ntripletsJGv)- ...
	     		     scl.K*sum((prefac(ajvolM_l1).*(CJ02+CJ13+CJ23)+Eps(ajvolM_l1)/(0.5*dt)).*ajvolS_l1)/linkL(l1);
	ntripletsJGv=ntripletsJGv+1;
	rowJGv(ntripletsJGv)=l1;
	colJGv(ntripletsJGv)=n2;
	valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+ ...
	     		     scl.K*sum((prefac(ajvolM_l1).*(CJ02+CJ13+CJ23)+Eps(ajvolM_l1)/(0.5*dt)).*ajvolS_l1)/linkL(l1);

	ntripletsJGH=ntripletsJGH+1;
	rowJGH(ntripletsJGH)=l1;
	colJGH(ntripletsJGH)=l1;
	valJGH(ntripletsJGH)=valJGH(ntripletsJGH)- ...
	     		     scl.K*sum((-prefac(ajvolM_l1).*(CJ02+CJ13+CJ23)-Eps(ajvolM_l1)/(0.5*dt)).*ajvolS_l1);
     end

JG_v=sparse2(rowJGv(1:ntripletsJGv),colJGv(1:ntripletsJGv),valJGv(1:ntripletsJGv),Nlink,Nnode);
JG_H=sparse2(rowJGH(1:ntripletsJGH),colJGH(1:ntripletsJGH),valJGH(1:ntripletsJGH),Nlink,Nlink);
clear rowJGv colJGv valJGv rowJGH colJGH valJGH;
JG_v = JG_v(eqnLinks,eqnNodes);
JG_H  = JG_H(eqnLinks,eqnLinks);
%savefilename='JG.mat';
%save(savefilename, 'JG_v','JG_H','-v7.3');
JG = [JG_v,JG_H];
clear JG_v JG_H;
    

display(['time for matrix collection,JG matrix:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jacob =sparse2([JF;JG]);

clear JF JG;

%%%0 preconditioned
%colind=colamd(Jacob);
%aspas = sparse2(Jacob(:,colind));
%[Lmatrix,Umatrix] = ilu(aspas,struct('type','ilutp','droptol',1e-4));	%most time consuming
%%%1. normal lu decomposition
[JacobLU.L, JacobLU.U, JacobLU.P, JacobLU.Q, JacobLU.R] = lu (Jacob) ;
%%%%2. klu decomposition from suitsparse
%[JacobLU,info] = klu (Jacob, struct ('ordering',0, 'scale', 1)) ;
%%%%3. cholmod from suitsparse
%[JacobLU,info] = ldlchol (Jacob);  
%fprintf ('ldlchol info : %d \n',info); 
%4 Factorize
%JacobLU = factorize(Jacob);
%%%5. normal lu decomposition without R or %%%6 LUsubs  
%[JacobLU.L, JacobLU.U, JacobLU.P, JacobLU.Q] = lu (Jacob) ;

%savefilename='JacobLU.mat';
%save(savefilename, 'JacobLU','-v7.3');

display(['time for matrix Jacob matrix combination:']);
toc;

toc(tStart);   
display(['  End building Jacob matrix.']);
