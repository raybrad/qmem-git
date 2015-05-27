function dJ_ci = dJ1dc1i(alpha,beta,type,vi)

switch type
    case 'n'
        dJ_ci = -vi*alpha*B(-beta);
    case 'p'
        dJ_ci = +vi*alpha*B(+beta);
    otherwise
        error('Invalid type of carrier');
    
end