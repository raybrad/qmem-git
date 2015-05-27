function  [lencr,rgncr] = semipnc(V,n,p)
%function  [currln,currlp,currrn,currrp] = semipnc(V,n,p)

global kx ky kz;
global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;
global sigma epsilon_in epsilon_sd;
global mun mup;
global scl;
global nodes links contacts;
global Nnode;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;

global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes dirNodes;
global isBndNodes isDirNodes dcVolDirNodes; %#ok<*NUSED>
global doping; % doping profile
global semiNodes eqnSemiNodes; %#ok<NUSED>
global dirSemiNodes isDirSemiNodes;

currln=0;
currlp=0;
currrn=0;
currrp=0;
totsareal= 0;
totsarear= 0;

for yi=qmy1:qmy2
   for zi=qmz1:qmz2

        n1=co2id([qmx1-1,yi,zi]);
        n2=co2id([qmx1,yi,zi]);;
        adlk = nodeLinks{n1}(1,:);
        adnd = nodeLinks{n1}(2,:);
        for j= 1:length(adlk)
            if adnd(j) == n2
               slink = adlk(j);
            end
        end
        v1  = V(n1);
        v2  = V(n2);
        ns1 = n(n1);
        ns2 = n(n2);
        ps1 = p(n1);
        ps2 = p(n2);

        slkvol  = linkVolumes{slink}(1,:);
        slkvols = linkVolumes{slink}(2,:);
        slkvolm = volumeM(slkvol);
        slinkl  = linkL(slink);
        alpha_n = mun/slinkl; alpha_p = mup/slinkl;

        for t = 1:length(slkvol)
            if slkvolm(t)==1
               beta = v2-v1;
               currln = currln + slkvols(t) * Jc(alpha_n,beta,'n',ns1,ns2);
               currlp = currlp + slkvols(t) * Jc(alpha_p,beta,'p',ps1,ps2);               
               totsareal= totsareal + slkvols(t);
            end
        end


        n1=co2id([qmx2,yi,zi]);
        n2=co2id([qmx2+1,yi,zi]);;
        adlk = nodeLinks{n1}(1,:);
        adnd = nodeLinks{n1}(2,:);
        for j= 1:length(adlk)
            if adnd(j) == n2
               slink = adlk(j);
            end
        end
        v1  = V(n1);
        v2  = V(n2);
        ns1 = n(n1);
        ns2 = n(n2);
        ps1 = p(n1);
        ps2 = p(n2);

        slkvol  = linkVolumes{slink}(1,:);
        slkvols = linkVolumes{slink}(2,:);
        slkvolm = volumeM(slkvol);
        slinkl  = linkL(slink);
        alpha_n = mun/slinkl; alpha_p = mup/slinkl;

        for t = 1:length(slkvol)
            if slkvolm(t)==1
               beta = v2-v1;
               currrn = currrn + slkvols(t) * Jc(alpha_n,beta,'n',ns1,ns2);
               currrp = currrp + slkvols(t) * Jc(alpha_p,beta,'p',ps1,ps2);
               totsarear= totsarear + slkvols(t);
            end
        end
   end
end

if (totsareal ~= 0)
    currln = currln/totsareal;
    currlp = currlp/totsareal;
end

if (totsarear ~= 0)
    currrn = currrn/totsarear;
    currrp = currrp/totsarear;
end

lencr = currln/(currln+currlp);
rgncr = currrn/(currrn+currrp);
