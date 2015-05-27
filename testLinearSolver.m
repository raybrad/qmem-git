%%%0 preconditioned
%colind=colamd(Jacob);
%aspas = sparse2(Jacob(:,colind));
%[Lmatrix,Umatrix] = ilu(aspas,struct('type','ilutp','droptol',1e-4));	%most time consuming
%%%1. normal lu decomposition
[JacobLU.L, JacobLU.U, JacobLU.P, JacobLU.Q, JacobLU.R] = lu (Jacob) ;
%%%%2. klu decomposition from suitsparse ..slow
%[JacobLU,info] = klu (Jacob, struct ('ordering',0, 'scale', 1)) ;
%%%%3. cholmod from suitsparse .. not working
%[JacobLU,info] = ldlchol (Jacob);  
%fprintf ('ldlchol info : %d \n',info); 
%4 Factorize  .. not as fast as 1
%JacobLU = factorize(Jacob);
%%%5. normal lu decomposition without R (..sim to 1) or %%%6 LUsubs( .. not working)
%[JacobLU.L, JacobLU.U, JacobLU.P, JacobLU.Q] = lu (Jacob) ;


%0
%[dX(colind),flag,relres,Itr]=bicgstab(Jacob(:,colind),rhs,1e-4,200,Lmatrix,Umatrix);
%display(['flag:',num2str(flag),' relres:',num2str(relres),' Itr:',num2str(Itr)]);
%1. normal lu decomposition
dX = JacobLU.Q * (JacobLU.U \ (JacobLU.L \ (JacobLU.P * (JacobLU.R \ rhs)))) ;
%2. klu decomposition from suitsparse ..slow
%dX = klu (JacobLU, '\', rhs) ;
%3. cholmod from suitsparse .. not working
%dX= ldlsolve(JacobLU,rhs);
%4 factorization  .. not as fast as 1
%dX = JacobLU\rhs ; 
%%%5. normal lu decomposition without R   almost similar to 1
%dX = JacobLU.Q * (JacobLU.U \ (JacobLU.L \ (JacobLU.P * rhs))) ;
%6 LUsubs,supposed to fasten the \ procedure  .. not working
%dX = LUsubs(JacobLU.L,JacobLU.U,JacobLU.P,JacobLU.Q,rhs);
