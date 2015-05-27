function [V,n,p,Vd,nd,pd] = staticSolution()

global isehpgr;


[Vs,ns,ps] = staticSolution;

  Vd = Vs;
  nd = ns;
  pd = ps;

if isehpgr == 1
  [Vs,ns,ps] = staticSolutionp(Vs,ns,ps);
end

  V=Vs;
  n=ns; 
  p=ps;
