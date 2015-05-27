function x = mylusolve(a,b,mtd)
switch mtd
    case 1
      x=a\b;
    case 2
      x=zeros(size(b));
      setup.droptol=1e-4;
      tol   =1e-4;
      maxit = 200;
      colind=colamd(a);
      aspas = sparse(a(:,colind));
%      [L,U] = ilu(aspas,struct('type','ilutp','droptol',1e-8));	%most time consuming
      [L,U] = ilu(aspas,struct('type','ilutp','droptol',1e-4));	%most time consuming
      x(colind)=bicgstab(a(:,colind),b,tol,maxit,L,U);
    case 3
%      path(['/home/lymeng/compecon/CEtools'],path);
%      path(['/home/lymeng/compecon/CEdemos'],path);
      colind=colamd(a);
      [L,U,P]=lu(a(:,colind));
      [rowind,i]=find(P');
      x=lusolve(L,U,b,rowind,colind);
    case 4
      x=zeros(size(b));
      setup.droptol=1e-4;
      tol   =1e-4;
      maxit = 200;
      aspas = sparse(a);
      [L,U] = ilu(aspas,setup);
      x=bicgstab(a,b,tol,maxit,L,U);
    otherwise
      error('undefined method for linear system');
end
