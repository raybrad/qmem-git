function [V,A,H,mJ_0,mJ_1,mJ_2,mJ_1p,mJ_2p,Efield_p,dtV,dtH] = tdsolverhs_revised(V_p,A_p,H_p,Js,mJ_0p,mJ_1p,mJ_2p,mJ_1pp,mJ_2pp,Efield_pp,dtV,dtH,dt)

display(['  Start TD simulation relaxation2 per step:']);

global omega_p gamma_p;
global lepsr_1 lomega_1 lgamma_1 lepsr_2 lomega_2 lgamma_2;
global links;
global Nnode Nlink;
global linkL; 
global dirNodes;
global metalLinks;
global EsurfLinks;
%global Jacob colind Lmatrix Umatrix;
global Fc11 Fc12 Fc21 Fc22 Fc30 Fn1matrix Fn2matrix Flkmatrix;
global Gc11 Gajlkmatrix11 Gc12 Gc20 Gn1matrix Gc31 Gajlk_n1matrix Gc32;
global Gn2matrix Gc41 Gajlk_n2matrix Gc42 Gc51 Gc52 Gc53 Gc54;
global JacobLU;
%Method 2
%Solving Gauss' law.
%Solving Current-continuity equations.
%Solving Maxwell-Ampere equations(already insert Lorentz gauge)
%Respecting the gauge condition is a side product of solving above set.

%%%%%%% initial guess  %%%%%%%%%%%%%
V = V_p + dtV *dt;
H = H_p + dtH *dt;
A = A_p + H_p *dt + dtH * dt^2/2;

mJ_0=mJ_0p;
mJ_1=mJ_1p;
mJ_2=mJ_2p;
Efield=Efield_pp;
Efield_p=Efield_pp;

dtV1 = dtV;
dtH1 = dtH;

eqnNodes = setdiff((1:Nnode)',dirNodes);%eqnNodes (nodes except from contact ,different from tdcalupdatec.m)   
					%dirNodes (contact part update by dtVp)

eqnLinks=setdiff((1:Nlink)',EsurfLinks);	%for E boundary condition

%%%%%%%% start Newton's iteration %%%%%%%%
updateTol = 1e-4;
droptol = 1e-6;
linSolveTol = 1e-6;
maxNewtonIt = 50;
maxLinIt = 200;
prefac=[1,0];

NeqnNodes = length(eqnNodes);
Nl = length(eqnLinks);

normRes = 1;
normUpdate = 1;
itNr = 0;
tStart=tic;

while itNr < maxNewtonIt && normUpdate > updateTol

tic;
   %update mJ,for plasmonic metal
   %drude mJ_(n+1)=(2-gamma_p*dt)/(2+gamma_p*dt)*mJ_n+ omega_p^2*dt/(2+gamma_p*dt)*(E_n+1 + E_n)
   %Lorentz mJ_(n+1)=alpha*mJ_n+beta*mJ_(n-1) + gamma * (E_n+1 -E_n-1)/2/dt 
   %or to save the storage, Lorentz mJ_(n+1)=alpha*mJ_n+beta*mJ_(n-1) + gamma* (E_n+1 -E_n)/dt 
   %alpha=(2-lomega*lomega*dt*dt)/(1+lgamma*dt);beta=(lgamma*dt-1)/(1+lgamma*dt);gamma=epsr*lomega*lomega*dt*dt/(1+lgamma*dt)
   if (~isempty(metalLinks))
   Efield(metalLinks)  = -(V(links(metalLinks,2))-V(links(metalLinks,1)))./linkL(metalLinks)-H(metalLinks);
   Efield_p(metalLinks)= -(V_p(links(metalLinks,2))-V_p(links(metalLinks,1)))./linkL(metalLinks)-H_p(metalLinks);
   
 %constant coefficient for plasmonic metal
   CJ01=(2-gamma_p*dt)/(2+gamma_p*dt);
   CJ02=dt*omega_p*omega_p/(2+gamma_p*dt);

   CJ11=(2-lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt);
   CJ12=(lgamma_1*dt-1)/(1+lgamma_1*dt);
   CJ13=(lepsr_1*lomega_1*lomega_1*dt*dt)/(1+lgamma_1*dt)/(2.0*dt);

   CJ21=(2-lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt);
   CJ22=(lgamma_2*dt-1)/(1+lgamma_2*dt);
   CJ23=(lepsr_2*lomega_2*lomega_2*dt*dt)/(1+lgamma_2*dt)/(2.0*dt);

   mJ_0(metalLinks)=CJ01.*mJ_0p(metalLinks)+CJ02.*(Efield(metalLinks)+Efield_p(metalLinks));

   mJ_1(metalLinks)=CJ11.*mJ_1p(metalLinks)+CJ12.*mJ_1pp(metalLinks)+CJ13.*(Efield(metalLinks)-Efield_pp(metalLinks));

   mJ_2(metalLinks)=CJ21.*mJ_2p(metalLinks)+CJ22.*mJ_2pp(metalLinks)+CJ23.*(Efield(metalLinks)-Efield_pp(metalLinks));
   end
display(['time for updating metalLinks:']);
toc;
tic;

 % construct the sparse matrix by (Row,Column,Value) 
rhs_F = zeros(Nnode,1);
FdtVmatrix=zeros(Nnode,6);
FdtHmatrix=zeros(Nnode,6);
FVmatrix  =zeros(Nnode,6);
FHmatrix  =zeros(Nnode,6);
FJmatrix  =zeros(Nnode,6);

FdtVmatrix(eqnNodes,:)= dtV(Fn2matrix(eqnNodes,:))-dtV(Fn1matrix(eqnNodes,:)); 
FdtHmatrix(eqnNodes,:)= dtH(Flkmatrix(eqnNodes,:));
FVmatrix(eqnNodes,:)  = V(Fn2matrix(eqnNodes,:))-V(Fn1matrix(eqnNodes,:));
FHmatrix(eqnNodes,:)  = H(Flkmatrix(eqnNodes,:));
FJmatrix(eqnNodes,:)  = (mJ_0(Flkmatrix(eqnNodes,:))+mJ_1(Flkmatrix(eqnNodes,:))+mJ_2(Flkmatrix(eqnNodes,:)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs_F=sum((Fc11.*FdtVmatrix+Fc12.*FdtHmatrix+Fc21.*FVmatrix+Fc22.*FHmatrix+Fc30.*FJmatrix),2);

display(['time for rhs_F matrix collection:']);
toc;
tic;


rhs_G =zeros(Nlink,1);
GAmatrix11=zeros(Nlink,12);
GAmatrix31=zeros(Nlink,6);
GdtVmatrix32=zeros(Nlink,1);
GAmatrix41=zeros(Nlink,6);
GdtVmatrix42=zeros(Nlink,1);
GdtVmatrix51=zeros(Nlink,1);

GAmatrix11(eqnLinks,:)=(Gajlkmatrix11(eqnLinks,:)~=0).*A(Gajlkmatrix11(eqnLinks,:)+(Gajlkmatrix11(eqnLinks,:)==0));

GAmatrix31(eqnLinks,:)=(Gajlk_n1matrix(eqnLinks,:)~=0).*A(Gajlk_n1matrix(eqnLinks,:)+(Gajlk_n1matrix(eqnLinks,:)==0));
GdtVmatrix32(eqnLinks)=dtV(Gn1matrix(eqnLinks));

GAmatrix41(eqnLinks,:)=(Gajlk_n2matrix(eqnLinks,:)~=0).*A(Gajlk_n2matrix(eqnLinks,:)+(Gajlk_n2matrix(eqnLinks,:)==0));
GdtVmatrix42(eqnLinks)=dtV(Gn2matrix(eqnLinks));
GdtVmatrix51(eqnLinks)=dtV(Gn2matrix(eqnLinks))-dtV(Gn1matrix(eqnLinks));

rhs_G = sum(Gc11.*GAmatrix11,2)+Gc12.*A+Gc20.*H+sum(Gc31.*GAmatrix31,2)+Gc32.*GdtVmatrix32+...
	sum(Gc41.*GAmatrix41,2)+Gc42.*GdtVmatrix42+Gc51.*GdtVmatrix51+Gc52.*dtH+Gc53.*(mJ_0+mJ_1+mJ_2)+Gc54.*Js;  
display(['time for rhs_G matrix collection:']);
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rhs_F = rhs_F(eqnNodes);
    rhs_G = rhs_G(eqnLinks);
      
    rhs = -sparse([rhs_F;rhs_G]);

    display(['time for matrix reconstruction:']);
    toc;
    tic;

    clear rhs_F rhs_G;
%     Jacob1 = sparse(Jacob);
%     rhs1 = sparse(rhs);
%     [L,U] = luinc(Jacob1,droptol);
%     [dX,flag,relres,iter] = gmres(Jacob1,rhs1,[],linSolveTol,maxLinIt,L,U);
%     if flag ~= 0, error('No convergence for linear solver'); end
%     gmresItNr(l) = iter(2);
%    invJacob = Jacob^(-1);
%    dX = invJacob*rhs;
%    dX = Jacob\rhs; newton raphson: x_n+1 = x_n - f(x_n)/f'(x_n)

   % dX = mylusovle(Jacob,rhs,1);
    dX=zeros(size(rhs));
    %[dX(colind),flag,relres,Itr]=bicgstab(Jacob(:,colind),rhs,1e-4,200,Lmatrix,Umatrix);
    %display(['flag:',num2str(flag),' relres:',num2str(relres),' Itr:',num2str(Itr)]);
    dX = JacobLU.Q * (JacobLU.U \ (JacobLU.L \ (JacobLU.P * (JacobLU.R \ rhs)))) ;

    display(['time for solving linear equation:']);
    toc;
    tic;

    dV = dX(1:NeqnNodes);
    dH = dX(NeqnNodes+1:NeqnNodes+Nl);
    
    normRes = norm(rhs);

    itNr = itNr+1;
    
    V(eqnNodes) = V(eqnNodes)+dV;
    H(eqnLinks) = H(eqnLinks)+dH;

    dtV  = 2 * (V - V_p)/dt - dtV1;
    dtH  = 2 * (H - H_p)/dt - dtH1;

    A(eqnLinks) = A_p(eqnLinks) + (H(eqnLinks) + H_p(eqnLinks))*dt/2;%+dtH(eqnLinks)*dt^2/2;
    display(['time for variables change']);
    toc;

    normUpdate = max([normRes]);%,norm(dV)/norm(V(eqnNodes)),norm(dH)/norm(H(eqnLinks))]);

    maxdV=max([norm(dV)/norm(V(eqnNodes))]);
    maxdH=max([norm(dH)/norm(H(eqnLinks))]);
    
    display([' maxdV: ',num2str(maxdV)]);
    display([' maxdH: ',num2str(maxdH)]);

    display(['  Iter:',num2str(itNr),'; ','normRes:',num2str(normRes),...
        '; ','normUpdate:',num2str(normUpdate),'; ','normdV:',num2str(norm(dV))]);

    maxmJ0=max([mJ_0]);
    maxmJ1=max([mJ_1]);
    maxmJ2=max([mJ_2]);
    display(['max mJ0:',num2str(maxmJ0)]);
    display(['max mJ1:',num2str(maxmJ1)]);
    display(['max mJ2:',num2str(maxmJ2)]);
    clear dX dV dH;
   if (itNr - maxNewtonIt == 0) && normUpdate > updateTol
	display(['reach maxNewtonIt and still not converge']);
   end 
    
end    %( end of while) 
toc(tStart);   
display(['  End TD simulation relaxation2 per step.']);
