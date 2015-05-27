function dJ_ci = dJdci(alpha,beta,type)

switch type
    case 'n'
        dJ_ci = -alpha*B(-beta);
    case 'p'
        dJ_ci = +alpha*B(+beta);
    otherwise
        error('Invalid type of carrier');
    
end