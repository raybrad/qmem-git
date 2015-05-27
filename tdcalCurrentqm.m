function [Jx,Jy,Jz] = calCurrent(nodes,V,n,p,A,H,dtVp,dtHp,iscd)

%%% calculate the six components of current at the given nodes %%%%%%% 
%%% Jx contains the currents on the two links at the left (minus ) and right (plus) sides of the node 

global sigma epsilon_in mun mup;
global volumeM linkL links linkS;
global nodeLinks linkVolumes nodeVolumes;
global isDirSemiNodes doping;
global epsilon_sd;
global epsilon_mt;
global epsilon_qm;
global isqmvolm;

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
    sign_n1 = sign(ajnd_n1-n1);
    Jlk = zeros(length(ajlk_n1),1);
    dV1 = deltaV(doping(n1));
    for i = 1:length(ajlk_n1)
%display([num2str(k),' ', num2str(i)]);
        n2 = ajnd_n1(i);
        lk = ajlk_n1(i);
        ajvol_lk = linkVolumes{lk}(1,:);
        ajvolS_lk = linkVolumes{lk}(2,:);
        ajvolM_lk = volumeM(ajvol_lk);
        for qmi = 1:length(ajvol_lk)
          if isqmvolm(ajvol_lk(qmi))
           ajvolM_lk(qmi)=4;
          end
        end
%ajvol_lk
        alpha_n = mun/linkL(lk); alpha_p = mup/linkL(lk);

        dtE1 = -(dtVp(n2)-dtVp(n1))/linkL(lk) - sign_n1(i)*dtHp(lk); 

        if isDirSemiNodes(n2)
            dV2 = deltaV(doping(n2));
        else
            dV2 = 0;
        end
        J_ = zeros(length(ajvol_lk),1);
        for j = 1:length(ajvol_lk)
            switch ajvolM_lk(j)
                case 1
                    beta = V(n2)-V(n1)+sign_n1(i)*H(lk)*linkL(lk);
                    J_(j) = iscd*(Jc(alpha_n,beta,'n',n(n1),n(n2))...
                        +Jc(alpha_p,beta,'p',p(n1),p(n2)))+ epsilon_sd*dtE1;
                case 2
                    J_(j) = iscd*sigma*((V(n1)-V(n2)+dV1-dV2)/linkL(lk) - sign_n1(i)*H(lk))+ epsilon_mt*dtE1;
                case 3
                    J_(j) = epsilon_in*dtE1;
                case 4
                    beta = V(n2)-V(n1)+sign_n1(i)*H(lk)*linkL(lk);
                    J_(j) = iscd*(Jc(alpha_n,beta,'n',n(n1),n(n2))...
                        +Jc(alpha_p,beta,'p',p(n1),p(n2)))+ epsilon_qm*dtE1;
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
