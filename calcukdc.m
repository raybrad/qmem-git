function ehpkdns = calcukd(V)

global sigma epsilon_in epsilon_sd;
global epsilon_mt;
global mun mup;
global scl;
global nodes links contacts;
global Nnode Nlink Nsurf;
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

global isehpgr;
global ehpgr;
global ehpkd ehpkf ehpkr;
global ehpda;
global ebpkt;
global fedstrgb;

 H= zeros(Nlink,1);
 ehpkdns = zeros(Nnode,1);
 abe=calcuaec(V,H);

 for n1 = semiNodes.'
     pb    = fedstrgb*abe(n1);
     kdcpb = 1 + pb + pb^2/3;
     ehpkdns(n1)=ehpkd*kdcpb;
 end






