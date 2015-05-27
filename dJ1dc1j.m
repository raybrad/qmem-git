function dJ_cj = dJ1dc1j(alpha,beta,type,vi)

switch type
    case 'n'
        dJ_cj = +vi*alpha*B(+beta);        
    case 'p'      
        dJ_cj = -vi*alpha*B(-beta);
    otherwise
        error('Invalid type of carrier');
    
end
