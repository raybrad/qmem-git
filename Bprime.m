function b = Bprime(beta)

if beta > 1e2 
    b = 0;
elseif abs(beta) < 1e-3 
    b  = -0.5;
else
    b = ((1-beta)*exp(beta)-1)/(exp(beta)-1)^2;
end
