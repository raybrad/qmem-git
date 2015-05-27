function [A,H] = tdstaticaqm(V,n,p,Ae,He)

display(['  Start TD calculation for staic qm A:']);

global sigma epsilon_in epsilon_sd;
global epsilon_mt;
global mun mup;
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

%%%%%%% initial guess  %%%%%%%%%%%%%
A = Ae;
H = zeros(Nlink,1);

dtV = zeros(Nnode,1);
dtH = zeros(Nlink,1);

eqnLinks = setdiff(QMlinks,bndLinks);

%%%%%%%% start Newton's iteration %%%%%%%%
newtonTol = 1e-15;
updateTol = 1e-10;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 10;
maxLinIt = 200;
Eps = [epsilon_sd,epsilon_mt,epsilon_in];
Sgm = [0,sigma,0];

Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
normRes_pre = 1e10;
itNr = 0;
tic;
    
while itNr < maxNewtonIt && normUpdate > updateTol

     JG_A = zeros(Nlink,Nlink); 
          
     rhs_G = zeros(Nlink,1);
     

    %%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
     for l1 = eqnLinks'
        n1 = links(l1,1);
        n2 = links(l1,2);
        n3 = links(l1,3);
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
            coefV_n1 = scl.K*(sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1))*linkS(l1)/linkL(l1);

            rhs_G(l1) = rhs_G(l1)-(-sign_n1.*coef_n1')*A(ajlk_n1);
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
            coefV_n2 = scl.K*(sum(Eps(ajvolM_n2).*ajvolV_n2)/nodeV(n2))*linkS(l1)/linkL(l1);


            rhs_G(l1) = rhs_G(l1)-sign_n2.*coef_n2'*A(ajlk_n2);
        end

        JG_A(l1,:) = CC-GD;

        %%%%%% source current %%%%%%%%%%%
        if n3==1
            rhs_G(l1) = rhs_G(l1) - scl.K * linkS(l1)*currdlink(l1,4);
        else
            rhs_G(l1) = rhs_G(l1) + 0;
        end
     end
    
    
    
    JG_A  = JG_A(eqnLinks,eqnLinks);
    
    rhs_G = rhs_G(eqnLinks);

    JG = JG_A;

    clear JG_A;
      
    Jacob = JG;
    rhs = -rhs_G;
    
    clear JG;
    clear rhs_G ;

    dX = Jacob\rhs;
    dA = dX(1:Nl);

    normRes = norm(rhs);
    normRes_pre = normRes;

    itNr = itNr+1;
    
    A(eqnLinks) = A(eqnLinks)+dA;


    currt=currentStd([1,round(ky/2),round(kz/2)],V,n,p,H,dtV,dtH,1)-currentStd([kx,round(ky/2),round(kz/2)],V,n,p,H,dtV,dtH,1)+currentStd([round(kx/2),1,round(kz/2)],V,n,p,H,dtV,dtH,2)-currentStd([round(kx/2),ky,round(kz/2)],V,n,p,H,dtV,dtH,2)+currentStd([round(kx/2),round(ky/2),1],V,n,p,H,dtV,dtH,3)-currentStd([round(kx/2),round(ky/2),kz],V,n,p,H,dtV,dtH,3);

    currb = currentStd([2,round(ky/2),round(kz/2)],V,n,p,H,dtV,dtH,1);

    normUpdate = max([normRes,norm(dA)/norm(A(eqnLinks)),abs(currt)/abs(currb)]);

%    normUpdate = max([norm(dV)/norm(V(eqnNodes)),norm(dA)/norm(A(eqnLinks)),norm(dB)/norm(B(eqnLinks))]);


    display(['  Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
        '; ','normUpdate:',num2str(normUpdate)]);

    clear dA;
    
end    
toc;   

display(['  End TD calculation for staic qm A.']);

