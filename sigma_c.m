function sigma = sigma_c(mu_c,beta,type,c0i,c0j)

%%%%%%%%%% equivalent small signal conductivity of semiconductor %%%%%%%%%
switch type
    case 'n'
        sigma = -mu_c*(Bprime(-beta)*c0i+Bprime(+beta)*c0j);
    case 'p'
        sigma = -mu_c*(Bprime(+beta)*c0i+Bprime(-beta)*c0j);
    otherwise
        error('Invalid type of carrier');
end