function [Jx,Jy,Jz] = calCurrent(nodes,Vs,ns,ps)

%%% calculate the six components of current at the given nodes %%%%%%% 
%%% Jx contains the currents on the two links at the left (minus ) and right (plus) sides of the node 

global sigma epsilon_in mun mup;
global volumeM linkL links linkS;
global nodeLinks linkVolumes nodeVolumes;
global isDirSemiNodes doping;

Nn = length(nodes);
Jx = zeros(Nn,2);
Jy = zeros(Nn,2);
Jz = zeros(Nn,2);
for k = 1:Nn
    n1 = nodes(k);
    ajlk_n1 = nodeLinks{n1}(1,:);
    ajnd_n1 = nodeLinks{n1}(2,:);
    ajvol_n1 = nodeVolumes{n1}(1,:);
    ajvolV_n1 = nodeVolumes{n1}(2,:);
    ajvolM_n1 = volumeM(ajvol_n1);
    Jlk = zeros(length(ajlk_n1),1);
    dV1 = deltaV(doping(n1));
    for i = 1:length(ajlk_n1)
%display([num2str(k),' ', num2str(i)]);
        n2 = ajnd_n1(i);
        lk = ajlk_n1(i);
        ajvol_lk = linkVolumes{lk}(1,:);
        ajvolS_lk = linkVolumes{lk}(2,:);
        ajvolM_lk = volumeM(ajvol_lk);
%ajvol_lk
        alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);
        if isDirSemiNodes(n2)
            dV2 = deltaV(doping(n2));
        else
            dV2 = 0;
        end
        J_ = zeros(length(ajvol_lk),1);
        for j = 1:length(ajvol_lk)
            switch ajvolM_lk(j)
                case 1
                    beta = Vs(n2)-Vs(n1);
                    J_(j) = Jc(alpha_n,beta,'n',ns(n1),ns(n2))...
                        +Jc(alpha_p,beta,'p',ps(n1),ps(n2));
                case 2
                    J_(j) = sigma*(Vs(n1)-Vs(n2)+dV1-dV2)/linkL(lk);
                case 3
                    J_(j) = 0;
                otherwise
                    error('undefined material');
            end
            
        end
        Jlk(i) = sum(J_.*ajvolS_lk.')/linkS(lk)^0;  
        linkDim = links(lk,3);
        linkID = 1*(n2 < n1)+2*(n2 > n1); 
        switch linkDim
            case 1
                Jx(k,linkID) = Jlk(i); 
            case 2
                Jy(k,linkID) = Jlk(i); 
            case 3
                Jz(k,linkID) = Jlk(i);   
        end    
    end
    
end
sum(Jx+Jy+Jz)
