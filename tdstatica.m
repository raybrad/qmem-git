function [A,H] = tdstatica(V)

display(['  Start TD calculation for static A:']);

global epsilon_in;
global epsilon_mt;
global scl;
global nodes links;
global Nnode Nlink;
global nodeLinks linkSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL volumeM;
global isBndNodes;
global bndLinks;

%%%%%%% initial guess  %%%%%%%%%%%%%
A = zeros(Nlink,1);
H = zeros(Nlink,1);

dtV = zeros(Nnode,1);
dtH = zeros(Nlink,1);

eqnLinks = setdiff((1:Nlink)',bndLinks);

%%%%%%%% start Newton's iteration %%%%%%%%
updateTol = 1e-6;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 10;
maxLinIt = 200;
Eps = [epsilon_mt,epsilon_in];

Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
normRes_pre = 1e10;
itNr = 0;
tic;
    
while itNr < maxNewtonIt && normUpdate > updateTol

     ntripletsJGA=(12+4*4+30)*Nlink;      
     rowJGA=zeros(ntripletsJGA,1);  
     colJGA=zeros(ntripletsJGA,1);  
     valJGA=zeros(ntripletsJGA,1);  
     ntripletsJGA=0;            
          
     rhs_G = zeros(Nlink,1);
     

    %%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
    %G = cur curl A - {grad div A+grad V}+{-J+par E/par t} =0
     for l1 = eqnLinks'                                      
        n1 = links(l1,1);
        n2 = links(l1,2);
        ajvol_l1 = linkVolumes{l1}(1,:);%volumeid
        ajvolS_l1 = linkVolumes{l1}(2,:);%associate surf area,not dual area, but only its own 1/4
        ajvolM_l1 = volumeM(ajvol_l1);
        CC = zeros(1,Nlink);
        GD = zeros(1,Nlink);
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
	    ntripletsJGA=ntripletsJGA+1;
	    rowJGA(ntripletsJGA)=l1;
	    colJGA(ntripletsJGA)=ajlk(m);
	    valJGA(ntripletsJGA)=CC(ajlk(m));
	    end
        end
        rhs_G(l1) = rhs_G(l1)+CC(l1)*A(l1);
	ntripletsJGA=ntripletsJGA+1;
	rowJGA(ntripletsJGA)=l1;
	colJGA(ntripletsJGA)=l1;
	valJGA(ntripletsJGA)=CC(l1);

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

            rhs_G(l1) = rhs_G(l1)-(-sign_n1.*coef_n1')*A(ajlk_n1);
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


            rhs_G(l1) = rhs_G(l1)-sign_n2.*coef_n2'*A(ajlk_n2);
        end

	    ajlk_n1n2=unique([ajlk_n1;ajlk_n2]);
            for m=1:length(ajlk_n1n2)
	    ntripletsJGA=ntripletsJGA+1;
	    rowJGA(ntripletsJGA)=l1;
	    colJGA(ntripletsJGA)=ajlk_n1n2(m);
	    valJGA(ntripletsJGA)=-GD(ajlk_n1n2(m));
	    end

        %%%%%% source current %%%%%%%%%%%
        E1 = -(V(n2)-V(n1))/linkL(l1)-H(l1);
        dtE1 = 0;
	dV1=0;
	dV2=0;
        for i = 1:length(ajvol_l1)
            switch ajvolM_l1(i)
                case 1
                    J =0.0 + epsilon_mt*dtE1;
                case 2
                    J = epsilon_in*dtE1;
                otherwise
                    error('undefined material');
            end
            rhs_G(l1) = rhs_G(l1)-scl.K*ajvolS_l1(i)*J;
        end
     end
    
    
    
    JG_A=sparse(rowJGA(1:ntripletsJGA),colJGA(1:ntripletsJGA),valJGA(1:ntripletsJGA),Nlink,Nlink);%ntripletsJGH,ntripletsJGH);
    JG_A  = JG_A(eqnLinks,eqnLinks);
    
    rhs_G = rhs_G(eqnLinks);

    JG = JG_A;

    clear JG_A;
      
    Jacob = JG;
    rhs = -rhs_G;
    
    clear JG;
    clear rhs_G ;

    dX = Jacob\rhs; %inv(Jacob)*rhs
    dA = dX(1:Nl);

    normRes = norm(rhs);
    normRes_pre = normRes;

    itNr = itNr+1;
    
    A(eqnLinks) = A(eqnLinks)+dA;


    normUpdate = max([normRes]);


    display(['  Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
        '; ','normUpdate:',num2str(normUpdate)]);

    clear dA dX;
    clear Jacob rhs;
    
end    
toc;   

display(['  End TD calculation for static A.']);

