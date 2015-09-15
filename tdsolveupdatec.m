function [mJ_0,mJ_1,mJ_2,Efield_p,dtV,dtH] = tdsolveupdatec(V_p,V_pp,A_p,H_p,H_pp,Js,mJ_0p,mJ_1p,mJ_2p,mJ_1pp,mJ_2pp,Efield_pp,dtVp,dtHp,dt)
%function [dtV,dtH] = tdsolveupdatec(V_p,A_p,H_p,Js,dtVp,dtHp,dt)

display(['  Start TD simulation relaxation2 per step:']);

global CJ01 CJ02 CJ11 CJ12 CJ13 CJ21 CJ22 CJ23;
global links;
global Nnode Nlink;
global linkL; 
global dirNodes;
global metalLinks;
global EsurfLinks;
global currdlink;
global sQMlinks;
global Fn1matrix;
global Fc;
global Gajlkmatrix11 Gajlk_n1matrix Gajlk_n2matrix Gn1matrix Gn2matrix Gn3matrix;
global Gc;
%Method 2
%Solving Gauss' law.
%Solving Current-continuity equations.
%Solving Maxwell-Ampere equations(already insert Lorentz gauge)
%Respecting the gauge condition is a side product of solving above set.

%%%%%%% initial guess  %%%%%%%%%%%%%
dtV    = dtVp;
dtH    = dtHp;

mJ_0=mJ_0p;
mJ_1=mJ_1p;
mJ_2=mJ_2p;
Efield=Efield_pp;
Efield_p=Efield_pp;


eqnNodes = setdiff((1:Nnode)',dirNodes);%eqnNodes (nodes except from contact ,different from tdcalupdatec.m)   
					%dirNodes (contact part update by dtVp)

eqnLinks=setdiff((1:Nlink)',EsurfLinks);	%for E boundary condition

%%%%%%%% start Newton's iteration %%%%%%%%
prefac=[1,0];

NeqnNodes = length(eqnNodes);
Nl = length(eqnLinks);

tStart=tic;



tic;
   %update mJ,for plasmonic metal
   %drude mJ_(n+1)=(2-gamma_p*dt)/(2+gamma_p*dt)*mJ_n+ omega_p^2*dt/(2+gamma_p*dt)*(E_n+1 + E_n)
   %Lorentz mJ_(n+1)=alpha*mJ_n+beta*mJ_(n-1) + gamma * (E_n+1 -E_n-1)/2/dt 
   %or to save the storage, Lorentz mJ_(n+1)=alpha*mJ_n+beta*mJ_(n-1) + gamma* (E_n+1 -E_n)/dt 
   %alpha=(2-lomega*lomega*dt*dt)/(1+lgamma*dt);beta=(lgamma*dt-1)/(1+lgamma*dt);gamma=epsr*lomega*lomega*dt*dt/(1+lgamma*dt)
   if (~isempty(metalLinks))
   Efield(metalLinks)  = -(V_p(links(metalLinks,2))-V_p(links(metalLinks,1)))./linkL(metalLinks)-H_p(metalLinks);
   Efield_p(metalLinks)= -(V_pp(links(metalLinks,2))-V_pp(links(metalLinks,1)))./linkL(metalLinks)-H_pp(metalLinks);
   
   mJ_0(metalLinks)=CJ01.*mJ_0p(metalLinks)+CJ02.*(Efield(metalLinks)+Efield_p(metalLinks));

   mJ_1(metalLinks)=CJ11.*mJ_1p(metalLinks)+CJ12.*mJ_1pp(metalLinks)+CJ13.*(Efield(metalLinks)-Efield_pp(metalLinks));

   mJ_2(metalLinks)=CJ21.*mJ_2p(metalLinks)+CJ22.*mJ_2pp(metalLinks)+CJ23.*(Efield(metalLinks)-Efield_pp(metalLinks));
   end
display(['time for updating metalLinks:']);
toc;
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rhsF
%rhs_F = zeros(Nnode,1);
%FdtVmatrix=zeros(Nnode,6);
%FdtHmatrix=zeros(Nnode,6);
%FVmatrix  =zeros(Nnode,6);
%FHmatrix  =zeros(Nnode,6);
%FJmatrix  =zeros(Nnode,6);

%FdtVmatrix(eqnNodes,:)= dtV(Fn2matrix(eqnNodes,:))-dtV(Fn1matrix(eqnNodes,:)); 
%FdtHmatrix(eqnNodes,:)= dtH(Flkmatrix(eqnNodes,:));
%FVmatrix(eqnNodes,:)  = V(Fn2matrix(eqnNodes,:))-V(Fn1matrix(eqnNodes,:));
%FHmatrix(eqnNodes,:)  = H(Flkmatrix(eqnNodes,:));
%FJmatrix(eqnNodes,:)  = (mJ_0(Flkmatrix(eqnNodes,:))+mJ_1(Flkmatrix(eqnNodes,:))+mJ_2(Flkmatrix(eqnNodes,:)));

%Fmatrix=[FdtVmatrix,FdtHmatrix,FVmatrix,FHmatrix,FJmatrix];
%rhs_F=sum((Fc.*Fmatrix),2);
%rhs_F=sum((Fc11.*FdtVmatrix+Fc12.*FdtHmatrix+Fc21.*FVmatrix+Fc22.*FHmatrix+Fc30.*FJmatrix),2);

FAmatrix  =zeros(Nnode,6);
FAmatrix(eqnNodes,:)=A_p(Fn1matrix(eqnNodes,:));
dtV=sum((Fc.*FAmatrix),2);
display(['time for getting dtV:']);
toc;
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rhsG

%rhs_G =zeros(Nlink,1);
GAmatrix11=zeros(Nlink,12);
GAmatrix31=zeros(Nlink,6);
GdtVmatrix32=zeros(Nlink,1);
GAmatrix41=zeros(Nlink,6);
GdtVmatrix42=zeros(Nlink,1);
GdtVmatrix51=zeros(Nlink,1);
Gqmcurrdlink=zeros(Nlink,6);

%GAmatrix11(eqnLinks,:)=(Gajlkmatrix11(eqnLinks,:)~=0).*A_p(Gajlkmatrix11(eqnLinks,:)+(Gajlkmatrix11(eqnLinks,:)==0));
%GAmatrix11(eqnLinks,:)=A_p(Gajlkmatrix11(eqnLinks,:)+(Gajlkmatrix11(eqnLinks,:)==0));
GAmatrix11(eqnLinks,:)=A_p(Gajlkmatrix11(eqnLinks,:));

%GAmatrix31(eqnLinks,:)=(Gajlk_n1matrix(eqnLinks,:)~=0).*A_p(Gajlk_n1matrix(eqnLinks,:)+(Gajlk_n1matrix(eqnLinks,:)==0));
%GAmatrix31(eqnLinks,:)=A_p(Gajlk_n1matrix(eqnLinks,:)+(Gajlk_n1matrix(eqnLinks,:)==0));
GAmatrix31(eqnLinks,:)=A_p(Gajlk_n1matrix(eqnLinks,:));
GdtVmatrix32(eqnLinks)=dtVp(Gn1matrix(eqnLinks));

%GAmatrix41(eqnLinks,:)=(Gajlk_n2matrix(eqnLinks,:)~=0).*A_p(Gajlk_n2matrix(eqnLinks,:)+(Gajlk_n2matrix(eqnLinks,:)==0));
%GAmatrix41(eqnLinks,:)=A_p(Gajlk_n2matrix(eqnLinks,:)+(Gajlk_n2matrix(eqnLinks,:)==0));
GAmatrix41(eqnLinks,:)=A_p(Gajlk_n2matrix(eqnLinks,:));
GdtVmatrix42(eqnLinks)=dtVp(Gn2matrix(eqnLinks));
GdtVmatrix51(eqnLinks)=dtVp(Gn2matrix(eqnLinks))-dtVp(Gn1matrix(eqnLinks));

for i=sQMlinks.'
	Gqmcurrdlink(i,Gn3matrix(i)+3)= currdlink(i,Gn3matrix(i)+3);	%
end

%Gmatrix=[GAmatrix11,A,H,GAmatrix31,GdtVmatrix32,GAmatrix41,GdtVmatrix42,GdtVmatrix51,dtH,(mJ_0+mJ_1+mJ_2),Js];
%rhs_G=sum((Gc.*Gmatrix),2);

%rhs_G = sum(Gc11.*GAmatrix11,2)+Gc12.*A+Gc20.*H+sum(Gc31.*GAmatrix31,2)+Gc32.*GdtVmatrix32+...
%	sum(Gc41.*GAmatrix41,2)+Gc42.*GdtVmatrix42+Gc51.*GdtVmatrix51+Gc52.*dtH+Gc53.*(mJ_0+mJ_1+mJ_2)+Gc54.*Js;  
          
Gmatrix=[GAmatrix11,A_p,H_p,GAmatrix31,GdtVmatrix32,GAmatrix41,GdtVmatrix42,GdtVmatrix51,(mJ_0+mJ_1+mJ_2),Js,Gqmcurrdlink];
%tmpGsum=zeros(Nlink,1);
%tmpGsum=sum((Gc.*Gmatrix),2);
%dtH(eqnLinks)= -tmpGsum(eqnLinks);
dtH=-sum((Gc.*Gmatrix),2);
%dtH =-Gc52LU.Q * (Gc52LU.U \ (Gc52LU.L \ (Gc52LU.P * (Gc52LU.R \ tmpGsum)))).*(Gc52(eqnLinks)~=0);
%dtH =-sum(Gc11.*GAmatrix11+Gc12.*A+Gc20.*H+sum(Gc31.*GAmatrix31,2)+Gc32.*GdtVmatrix32+sum(Gc41.*GAmatrix41,2)+Gc42.*GdtVmatrix42+Gc51.*GdtVmatrix51+Gc53.*(mJ_0+mJ_1+mJ_2)+Gc54.*Js,2)/Gc52;  
display(['time for getting dtH:']);
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    
toc(tStart);   
display(['  End updateing dtV dtH.']);
