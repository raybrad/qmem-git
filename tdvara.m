function [VT,AT,HT] = tdvara(nsteps)

display(['Start TD variables analyzation']);

global scl;
global Nnode Nlink;
global savefile;


%%%%%%% initial guess  %%%%%%%%%%%%%
VT = zeros(Nnode,nsteps);
AT = zeros(Nlink,nsteps);
HT = zeros(Nlink,nsteps);

ntp  = 0;
tic;

while ntp < nsteps
   
    
 filename = [savefile,num2str(ntp),'.mat'];
 load(filename);
 VT(:,ntp+1)   = V(:,1);
 AT(:,ntp+1)   = A(:,1);
 HT(:,ntp+1)   = H(:,1);


    ntp = ntp + 1;
    
end    
toc;   

display(['End TD variables analyzation']);

