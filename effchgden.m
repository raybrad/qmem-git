function effcd = effchgden(H)

global epsilon_in;
global epsilon_mt;
global epsilon_qm;
global Nnode;
global nodeLinks linkVolumes nodeVolumes;
global volumeM;
global metalNodes;
global isqmvolm;

effcd = zeros(Nnode,1); 

NeqnmetalNodes=length(metalNodes);

Eps = [epsilon_mt,epsilon_in,epsilon_qm];

for k=1:NeqnmetalNodes
    n1=metalNodes(k);
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
	%%
       %for qmi = 1:length(ajvol_lk)
       %    if isqmvolm(ajvol_lk(qmi))
       %       ajvolM_lk(qmi)=3;
       %    end
       %end
        effcd(n1) = effcd(n1)+sign_n1(i)* H(lk)*sum(Eps(ajvolM_lk).*ajvolS_lk);
    end
%    semiV = sum(ajvolV_n1(ajvolM_n1 == 1));
%    effcd(n1) = effcd(n1)/semiV;	%average over volume
     metalV = sum(ajvolV_n1(ajvolM_n1 == 2));	%metalVolume
     if ~(metalV==0)
     effcd(n1) = effcd(n1)/metalV;
     end
end
