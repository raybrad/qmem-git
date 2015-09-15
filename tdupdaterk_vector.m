function [Vn,An,Hn,mJ_0n,mJ_1n,mJ_2n,Efield_nn,dtVn,dtHn,Js] = tdupdaterk_vector(V,V_p,A,H,H_p,mJ_0,mJ_1,mJ_2,mJ_1p,mJ_2p,Efield_p,dtVp,dtHp,dt,ntp)
%k1=f(tn,yn)
%k2=f(tn+h/2,yn+h/2*k1)
%k3=f(tn+h/2,yn+h/2*k2)
%k4=f(tn+h,  yn+h*k3  )
%tn+1=tn+h
%yn+1=yn+h/6*(k1+2k2+2k3+h4)

global epdf epdf2;
global extfelc;
global dirNodes dcVolDirNodes acVolDirNodes;
global nsteps;
global irkod;
global lightsource tlas tzero;
global scl;
global Nlink JLinks J_amp;

%%
tic;
Js=zeros(Nlink,1);
timec = (ntp-1)*dt;
timee = timec + dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emv  = V*scl.Vt;
emdt = dt*scl.tao*1e15;
currtime = timee*scl.tao*1e15;
printVoltage(emv,currtime,emdt);
%tn+h
switch extfelc
   case 1
     dtVp(dirNodes) = acVolDirNodes*exp(- timee/epdf)/epdf;
   case 2
     dtVp(dirNodes) = pi*acVolDirNodes*sin(2*pi*timee/epdf)/epdf;
   case 3
     dtVp(dirNodes) = 0.0;
   case 4
     dtVp(dirNodes) = -acVolDirNodes*exp((timee-nsteps*dt)/epdf)/epdf;
   otherwise
     error('undefined bias profile');
end
% a current density sheet is added to simulate the light source
%tn+h
switch lightsource
 case 1
 Js(JLinks)=J_amp*sin(2*pi*timee/epdf2);
 case 2
 Js(JLinks)=J_amp*sin(2*pi*timee/epdf2)*exp(-((timee-tzero)/tlas)^2.0);   %omega= 2*pi/1 ,E = hbar *omega =0.6582*2*3.14/1=4.1356eV ~ 300nm
 case 3
 if (timee<=7*tlas) 
 Js(JLinks)=J_amp*exp(-((timee-tzero)/tlas)^2.0);   %omega= 2*pi/1 ,E = hbar *omega =0.6582*2*3.14/1=4.1356eV ~ 300nm
 else
 Js(JLinks)=0.0;
 end
 otherwise
 error('undefined light source');
end

%[dtV4,dtH4] = tdsolveupdatec(Ve,Ae,He,Js,dtVp,dtHp,dt);
[mJ_0n,mJ_1n,mJ_2n,Efield_nn,dtVn,dtHn] = tdsolveupdatec(V,V_p,A,H,H_p,Js,mJ_0,mJ_1,mJ_2,mJ_1p,mJ_2p,Efield_p,dtVp,dtHp,dt);

Vn = V + dt* dtVn;
Hn = H + dt* dtHn;
An = A + dt*   Hn+dt^2*dtHn/2;
%or in another way
%(Vn+1-Vn)/dt=sum((Fc11.*FAmatrix),2);
toc;
display(['  End time simulation per step.']);
