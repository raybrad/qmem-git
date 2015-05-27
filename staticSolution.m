function [Vs] = staticSolution()

%%% static solution of Poisson and drift-diffusion (DD) equations 
%%% use Newton's method
%%% only electron density n is used (p is all zero at this moment)

display('Start static solution');

global sigma epsilon_in epsilon_mt;
global mun mup;
global scl;
global nodes links contacts;
global Nnode;
global nodeLinks linkSurfs surfNodes surfLinks volumeNodes volumeLinks...
    volumeSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL nodeM linkM linkCenter surfCenter volumeM;
global bndNodes edgeNodes dirNodes;
global isBndNodes isDirNodes dcVolDirNodes;
global doping; % doping profile
global kx ky kz;
global bndLinks;
global metalNodes;

global eqnNodes;


%%%%%%% generate initial guess (equilibrium solution) %%%%%%%%%%%%%
Vs = zeros(Nnode,1);
Vs(dirNodes) = 1*dcVolDirNodes;

eqnNodes = setdiff((1:Nnode)',dirNodes);


%%% numbers of Poisson and DD equations
Nn = length(eqnNodes);

%%%%%%%% start Newton's iteration %%%%%%%%
updateTol = 1e-6;
maxNewtonIt = 10;
Eps = [epsilon_mt,epsilon_in];
Sgm = [sigma,0];
normRes = 1;
normUpdate = 1;
itNr = 0;

tStart=tic;
while itNr < maxNewtonIt && normUpdate > updateTol
    tic;

     rhs_F = zeros(Nnode,1); 

     ntripletsJFv=48*Nnode;
     rowJFv=zeros(ntripletsJFv,1);
     colJFv=zeros(ntripletsJFv,1);
     valJFv=zeros(ntripletsJFv,1);
     ntripletsJFv=0;

   % Gauss's law  %integral(epsilon*laplace V-rho)=0
   % current-continuity div*J=0
    for n1 = eqnNodes.'
        ajlk_n1 = nodeLinks{n1}(1,:);
        ajnd_n1 = nodeLinks{n1}(2,:);
        ajvol_n1 = nodeVolumes{n1}(1,:);
        ajvolV_n1 = nodeVolumes{n1}(2,:);
        ajvolM_n1 = volumeM(ajvol_n1);
            for i = 1:length(ajlk_n1)
                n2 = ajnd_n1(i);
                lk = ajlk_n1(i);

                    dV2 = 0;

                ajvol_lk = linkVolumes{lk}(1,:);
                ajvolS_lk = linkVolumes{lk}(2,:);
                ajvolM_lk = volumeM(ajvol_lk);

                if any(ajvolM_n1 == 1)
                 tcoem = 0;
                else 
                 tcoem = 1;
                end

                coefV = sum((tcoem*Eps(ajvolM_lk)+Sgm(ajvolM_lk)).*ajvolS_lk)/linkL(lk);
		    ntripletsJFv=ntripletsJFv+1;
	  	    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n1;
		    valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+coefV;
	 	    ntripletsJFv=ntripletsJFv+1;
		    rowJFv(ntripletsJFv)=n1;
		    colJFv(ntripletsJFv)=n2;
		    valJFv(ntripletsJFv)=-coefV;
                dV2_vec = [1,1/2];
                rhs_F(n1) = rhs_F(n1)+sum((tcoem*Eps(ajvolM_lk)+Sgm(ajvolM_lk)).*...
                        (-(Vs(n2)-Vs(n1)+dV2_vec(ajvolM_lk).*dV2)).*ajvolS_lk)/linkL(lk);
            end
        end

JF_v=sparse(rowJFv(1:ntripletsJFv),colJFv(1:ntripletsJFv),valJFv(1:ntripletsJFv),Nnode,Nnode);
JF_v = JF_v(eqnNodes,eqnNodes);
display(['time for matrix collection,F matrix:']);
toc;
tic;

rhs_F = rhs_F(eqnNodes);
Jacob = [JF_v];
rhs = -[rhs_F];

    clear JF_v;
    clear rhs_F;

    dX = mylusovle(Jacob,rhs,1);
    display(['time for solving linear equation:']);
    toc;
    dV = dX(1:Nn);
    normRes = norm(rhs);
    
    itNr = itNr+1;
    
    Vs(eqnNodes) = Vs(eqnNodes)+dV;

    normUpdate = max([normRes]);

	display(['Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
	        '; ','normUpdate:',num2str(normUpdate),'; ','normdV:',num2str(norm(dV))]);

	clear dX dV;
    	clear Jacob rhs;
end  

toc(tStart);   

display('End static solution');
