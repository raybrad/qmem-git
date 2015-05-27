function b = B(beta)

if abs(beta) < 1e-3
    b = 1;
elseif abs(beta) > 1e2
    b = 0;
else
    b = beta/(exp(beta)-1);
end