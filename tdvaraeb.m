function [VT,nT,pT,AT,HT,ET,BT,ANT] = tdvara(nsteps)

display(['Start TD variables(including E & B) analyzation']);

global scl;
global Nnode Nlink;
global savefile;


%%%%%%% initial guess  %%%%%%%%%%%%%
VT = zeros(Nnode,nsteps);
nT = zeros(Nnode,nsteps);
pT = zeros(Nnode,nsteps);
AT = zeros(Nlink,nsteps);
HT = zeros(Nlink,nsteps);
ET = zeros(Nnode,nsteps);
BT = zeros(Nnode,nsteps);
ANT= zeros(Nnode,nsteps);

ntp  = 0;
tic;

while ntp < nsteps
   
    
 filename = [savefile,num2str(ntp),'.mat'];
 load(filename);
 VT(:,ntp+1)   = V(:,1);
 nT(:,ntp+1)   = n(:,1);
 pT(:,ntp+1)   = p(:,1);
 AT(:,ntp+1)   = A(:,1);
 HT(:,ntp+1)   = H(:,1);

 [E,B,AN] = tdcalcueb(V,n,p,A,H);

 ET(:,ntp+1)   = E(:,4);
 BT(:,ntp+1)   = B(:,4);
 ANT(:,ntp+1)  = AN(:,4);


    ntp = ntp + 1;
    
end    
toc;   

display(['End TD variables(including E & B) analyzation']);

