function dJ_cj = dJdcj(alpha,beta,type)

switch type
    case 'n'
        dJ_cj = +alpha*B(+beta);        
    case 'p'      
        dJ_cj = -alpha*B(-beta);
    otherwise
        error('Invalid type of carrier');
    
end