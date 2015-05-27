function dJn_v = dJndV(alpha_n,beta,n1,n2)


dJn_v = -alpha_n*Bprime(-beta)*n1-alpha_n*Bprime(+beta)*n2;


