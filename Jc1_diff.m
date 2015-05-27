function J = Jc1_diff(alpha,beta,type,c1i,c1j,vi,vj)

%%%%%%% small signal semiconductor current (diffusion part) %%%%%%%%%%%%%
switch type
    case 'n'
        J = -vi*alpha*B(-beta)*c1i+vj*alpha*B(+beta)*c1j;
    case 'p'
        J = +vi*alpha*B(+beta)*c1i-vj*alpha*B(-beta)*c1j;
    otherwise
        error('Invalid type of carrier');
end