function [V,n,p,A,H] = tdsimulSolution(dt,nsteps)

display(['Start TD dynamic solution']);

global nedrelax;
global irkod;
global savefile;
global Nlink;
global Nnode;
global extfelc;
global dirNodes dcVolDirNodes acVolDirNodes;
global epdf;

global QMnodes;
global sQMnodes;
global isSQMlinks;
global QMlinks;
global isqmnodes;

%%%%%%% initial guess  %%%%%%%%%%%%%
%V = Vs;
%n = ns;
%p = ps;
A = zeros(Nlink,1);
H = zeros(Nlink,1);

%[A,H] = tdstatica(V,n,p);

initials = 'tdinis.mat';

load(initials);

savefilename = [savefile,num2str(0),'.mat'];
save(savefilename, 'V', 'n', 'p');

if nedrelax==2
  fp = fopen('tdqmcurrd.dat','w');
  fclose(fp);
end

ntp  = 1; 
tic;

while ntp < nsteps + 1
   
 display(['Time step:',num2str(ntp)]);

 [Vu,nu,pu] = statupdatenmc(V,n,p,A,H,0,ntp);
 
 V = Vu;
 n = nu;
 p = pu;
    
 savefilename = [savefile,num2str(ntp),'.mat'];
 save(savefilename, 'V', 'n', 'p');

    ntp = ntp + 1;
    
end    
toc;   

display(['End TD dynamic solution']);
