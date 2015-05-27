function effcd = effchgden(H)

global sigma epsilon_in epsilon_sd;
global epsilon_mt;
global epsilon_qm;
global mun mup;
global scl;
global nodes links contacts;
global Nnode Nlink;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolS nodeVolV;
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

%%%%%%% initial guess  %%%%%%%%%%%%%
effcd = zeros(Nnode,1); 

%eqnSemiNodes = semiNodes;

%eqnNodes = setdiff((1:Nnode)',[dirNodes]);
%NeqnNodes = length(eqnNodes);
NeqnmetalNodes=length(metalNodes);
%%%%%%%% start Newton's iteration %%%%%%%%
newtonTol = 1e-15;
updateTol = 1e-10;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 10;
maxLinIt = 200;
Eps = [epsilon_sd,epsilon_mt,epsilon_in,epsilon_qm];
Sgm = [0,sigma,0];

%for n1 = eqnSemiNodes.'
for k=1:NeqnmetalNodes
    n1=metalNodes(k);
    ajlk_n1 = nodeLinks{n1}(1,:);
    ajnd_n1 = nodeLinks{n1}(2,:);
    ajvol_n1 = find(nodeVolV(n1,:));
    ajvolV_n1 = nodeVolV(n1,ajvol_n1);
    ajvolM_n1 = volumeM(ajvol_n1);
    sign_n1 = sign(ajnd_n1-n1);

    for i = 1:length(ajlk_n1)
        n2 = ajnd_n1(i);
        lk = ajlk_n1(i);
        ajvolS_lk = linkVolS(lk,ajvol_lk);
        ajvolM_lk = volumeM(ajvol_lk);
        ajvolM_lk = volumeM(ajvol_lk);
        for qmi = 1:length(ajvol_lk)
            if isqmvolm(ajvol_lk(qmi))
               ajvolM_lk(qmi)=4;
            end
        end
        effcd(n1) = effcd(n1)+sign_n1(i)* H(lk)*sum(Eps(ajvolM_lk).*ajvolS_lk);
    end
%    semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
%    effcd(n1) = effcd(n1)/semiV;	%average over volume
     metalV = sum(ajvolV_n1(ajvolM_n1 == 2));	%metalVolume
     if ~(metalV==0)
     effcd(n1) = effcd(n1)/metalV;
     end
end
