function [dtV,dtH] = tdcalupdatec(V,A,H,Js,Yp,Zp)

global epsilon_in;
global epsilon_mt;
global omega_p gamma_p;
global lepsr_1 lomega_1 lgamma_1 lepsr_2 lomega_2 lgamma_2;
global scl;
global nodes links;
global Nnode Nlink;
global nodeLinks linkSurfs linkVolumes nodeVolumes;
global nodeV linkL linkS dlinkL volumeM;
global bndNodes dirNodes;
global isBndNodes;
global bndLinks;
global metalLinks;

%%%%%%% initial guess  %%%%%%%%%%%%%
Y    = Yp;
Z    = Zp;


nobndNodes = setdiff((1:Nnode)',bndNodes);

eqnNodes  = setdiff(bndNodes,dirNodes);
eqnNodes1 = setdiff(nobndNodes,dirNodes); 

eqnLinks = setdiff((1:Nlink)',bndLinks);

%%%%%%%% start Newton's iteration %%%%%%%%
updateTol = 1e-6;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 10;
maxLinIt = 200;
Eps = [epsilon_mt,epsilon_in];
prefac=[1,0];

NeqnNodes = length(eqnNodes);
Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
normRes_pre = 1e10;
itNr = 0;

tStart=tic;
%%%% Lorezentz gauge%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% par               par     1                      % 
%---- V= div*A  --> ----V+ ---- integral0(A)dS=0   % 
%par t              par t   Vol                    %
%get Y (dV/dt)                                     %
%%%%%%% initial guess  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for n1 = eqnNodes1.'
    ajlk_n1 = nodeLinks{n1}(1,:);
    ajnd_n1 = nodeLinks{n1}(2,:);
    ajvol_n1 = nodeVolumes{n1}(1,:);
    ajvolV_n1 = nodeVolumes{n1}(2,:);
    ajvolM_n1 = volumeM(ajvol_n1);
    sign_n1 = sign(ajnd_n1-n1);
    coef_n1 = linkS(ajlk_n1)./nodeV(n1);
    Y(n1)  = -sign_n1.*coef_n1'*A(ajlk_n1);
    epsnode = sum(Eps(ajvolM_n1).*ajvolV_n1)/nodeV(n1);
    Y(n1) = Y(n1) / (scl.K * epsnode);
end
%
    
while itNr < maxNewtonIt && normUpdate > updateTol

     tic;
     rhs_F = zeros(Nnode,1); 
     rhs_G = zeros(Nlink,1);

     ntripletsJFv=12*Nnode;
     rowJFv=zeros(ntripletsJFv,1);
     colJFv=zeros(ntripletsJFv,1);
     valJFv=zeros(ntripletsJFv,1);
     ntripletsJFv=0;

     ntripletsJFH=6*Nnode;
     rowJFH=zeros(ntripletsJFH,1);
     colJFH=zeros(ntripletsJFH,1);
     valJFH=zeros(ntripletsJFH,1);
     ntripletsJFH=0;

     ntripletsJGv=24*Nlink;     
     rowJGv=zeros(ntripletsJGv,1);  
     colJGv=zeros(ntripletsJGv,1);
     valJGv=zeros(ntripletsJGv,1);  
     ntripletsJGv=0;                

     ntripletsJGH=(12+4*4+30)*Nlink;      
     rowJGH=zeros(ntripletsJGH,1);  
     colJGH=zeros(ntripletsJGH,1);  
     valJGH=zeros(ntripletsJGH,1);  
     ntripletsJGH=0;                
    
     % F= div J_tot = 0 (for node, no A) J_tot =J_f +epsilon part E/par t    
     % current continuity equation + gauss' law: div J_f+ d (div D)/dt = 0
     % => J_f* S+ eps *dE/dt *S =0
     % sigma*E + par E/par t (metal)                                 
     % -miu_nij/h_ij*B[-(Vj-Vi+PI_ij*h_ij)]*n_ij+miu_nij/h_ij*B[(Vj-Vi+PI_ij*h_ij)]*n_ij
     for k=1:NeqnNodes
        n1 = eqnNodes(k);
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

            dtE1 = -(Y(n2)-Y(n1))/linkL(lk)-sign_n1(i)*Z(lk);		%? how about W
		ntripletsJFv=ntripletsJFv+1;
		rowJFv(ntripletsJFv)=n1;
		colJFv(ntripletsJFv)=n1;
		valJFv(ntripletsJFv)=valJFv(ntripletsJFv)+sum(Eps(ajvolM_lk).*ajvolS_lk)/linkL(lk);

		ntripletsJFv=ntripletsJFv+1;
		rowJFv(ntripletsJFv)=n1;
		colJFv(ntripletsJFv)=n2;
		valJFv(ntripletsJFv)=valJFv(ntripletsJFv)-(isBndNodes(n2))*sum(Eps(ajvolM_lk).*ajvolS_lk)/linkL(lk);

		ntripletsJFH=ntripletsJFH+1;
		rowJFH(ntripletsJFH)=n1;
		colJFH(ntripletsJFH)=lk;
		valJFH(ntripletsJFH)=valJFH(ntripletsJFH)+(-sign_n1(i))*sum(Eps(ajvolM_lk).*ajvolS_lk);
                
		rhs_F(n1) = rhs_F(n1)+sum(Eps(ajvolM_lk)*dtE1.*ajvolS_lk);

        end
     end

     display(['time for matrix collection,F matrix:']);
     toc;
     tic;
    %%%%%%%%%%%%% Build Jacobian and rhs of G (Ampere's law) %%%%%%%%%%%%%
    % G = cur curl A - {J_f+par E/par t} =0 (for link A,PI)   
    %(par PI /par t)_n+1 = (par PI/par t)_n -G/{par G/par(par PI/par t)}   to get the (par PI/par t)
     for k=1:Nl
        l1 = eqnLinks(k);
        n1 = links(l1,1);
        n2 = links(l1,2);
        ajvol_l1 = linkVolumes{l1}(1,:);%volumeid
        ajvolS_l1 = linkVolumes{l1}(2,:);%associate surf area,not dual area, but only its own 1/4
        ajvolM_l1 = volumeM(ajvol_l1);
        CC = zeros(1,Nlink);

        %%%%%%%%%% curl curl of A %%%%%%%%%%
        ajsf_l1 = linkSurfs{l1};
        for i = 1:size(ajsf_l1,2)
            ajlk = ajsf_l1(2:4,i);
            S = linkL(l1)*linkL(ajlk(1));
            dL = dlinkL(ajsf_l1(1,i));
            ajlkL = linkL(ajlk)';
            if n2 == links(ajlk(1),1)
                ajlkSign = [1,-1,-1];
            elseif n2 == links(ajlk(1),2)
                ajlkSign = [-1,-1,1];
            else
                error('Incorrect link arrangement');
            end
            CC(ajlk) = dL/S*(ajlkSign.*ajlkL);
            CC(l1) = CC(l1)+dL/S*linkL(l1);
            rhs_G(l1) = rhs_G(l1)+CC(ajlk)*A(ajlk);
        end
        rhs_G(l1) = rhs_G(l1)+CC(l1)*A(l1);

        %%%%%% source current %%%%%%%%%%%
        dtE1 = -(Y(n2)-Y(n1))/linkL(l1) - Z(l1);

        ntripletsJGv=ntripletsJGv+1;
        rowJGv(ntripletsJGv)=l1;
	colJGv(ntripletsJGv)=n1;
	valJGv(ntripletsJGv)=valJGv(ntripletsJGv)-(isBndNodes(n1))*scl.K*sum(Eps(ajvolM_l1).*ajvolS_l1)/linkL(l1);

        ntripletsJGv=ntripletsJGv+1;
        rowJGv(ntripletsJGv)=l1;
	colJGv(ntripletsJGv)=n2;
	valJGv(ntripletsJGv)=valJGv(ntripletsJGv)+(isBndNodes(n2))*scl.K*sum(Eps(ajvolM_l1).*ajvolS_l1)/linkL(l1);


        ntripletsJGH=ntripletsJGH+1;
 	rowJGH(ntripletsJGH)=l1;
	colJGH(ntripletsJGH)=l1;
	valJGH(ntripletsJGH)=valJGH(ntripletsJGH)-scl.K*sum(-Eps(ajvolM_l1).*ajvolS_l1);

        rhs_G(l1) = rhs_G(l1)-scl.K*sum((Eps(ajvolM_l1)*dtE1+prefac(ajvolM_l1)*Js(l1)).*ajvolS_l1);
     end
display(['time for matrix collection,G matrix:']);
toc;
tic;
    
    JF_v=sparse2(rowJFv(1:ntripletsJFv),colJFv(1:ntripletsJFv),valJFv(1:ntripletsJFv),Nnode,Nnode);
    JF_H=sparse2(rowJFH(1:ntripletsJFH),colJFH(1:ntripletsJFH),valJFH(1:ntripletsJFH),Nnode,Nlink);
    clear rowJFv colJFv valJFv rowJFH colJFH valJFH;
    JF_v = JF_v(eqnNodes,eqnNodes);
    JF_H  = JF_H(eqnNodes,eqnLinks);
    
    JG_v=sparse2(rowJGv(1:ntripletsJGv),colJGv(1:ntripletsJGv),valJGv(1:ntripletsJGv),Nlink,Nnode);
    JG_H=sparse2(rowJGH(1:ntripletsJGH),colJGH(1:ntripletsJGH),valJGH(1:ntripletsJGH),Nlink,Nlink);
    clear rowJGv colJGv valJGv rowJGH colJGH valJGH;
    JG_v = JG_v(eqnLinks,eqnNodes);
    JG_H  = JG_H(eqnLinks,eqnLinks);
    
    rhs_F = rhs_F(eqnNodes);
    rhs_G = rhs_G(eqnLinks);

    JF = [JF_v,JF_H];
    JG = [JG_v,JG_H];

    clear JF_v JF_H JG_v JG_H;

    Jacob = [JF;JG];
    rhs = -[rhs_F;rhs_G];
    display(['time for matrix reconstruction:']);
    toc;
    tic;
  
    clear JF JG rhs_F rhs_G;
    
    dX = Jacob\rhs;
    dY = dX(1:NeqnNodes);
    dZ = dX(NeqnNodes+1:NeqnNodes+Nl);

    display(['time for solving linear equation:']);
    toc;

    normRes = norm(rhs);
    normRes_pre = normRes;

    itNr = itNr+1;
    
    Y(eqnNodes) = Y(eqnNodes)+dY;   %dtV
    Z(eqnLinks) = Z(eqnLinks)+dZ;   %dtPI
    maxdY=max([dY]);
    maxdZ=max([dZ]);
    display(['max dY:',num2str(maxdY)]);
    display(['max dZ:',num2str(maxdZ)]);

    normUpdate = max([normRes]);
    display(['  Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),'; ','normUpdate:',num2str(normUpdate)]);
    clear dX dY dZ;    
end    
%R
dtV = Y;
dtH = Z;

toc(tStart);   

display(['  End tdcalupdatec .']);
