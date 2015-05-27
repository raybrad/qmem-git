function dJp_v = dJpdV(alpha_p,beta,p1,p2)


dJp_v = -alpha_p*Bprime(beta)*p1-alpha_p*Bprime(-beta)*p2;
