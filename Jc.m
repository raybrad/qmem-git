function J = Jc(alpha,beta,type,c1,c2)

switch type
    case 'n'
        J = -alpha*B(-beta)*c1+alpha*B(+beta)*c2;
    case 'p'
        J = +alpha*B(+beta)*c1-alpha*B(-beta)*c2;
    otherwise
        error('Not valid type of carrier');
end