function dV = deltaV(N)

% contact potential between metal/semiconductor interface (ohmic contact)

global ni;

dV = zeros(length(N),1);
idp = N > 0;
idm = N < 0;
dV(idp) = -log(N(idp)./(2*ni).*(1+sqrt(1+4*ni^2./N(idp).^2)));
dV(idm) = log(-N(idm)./(2*ni).*(1+sqrt(1+4*ni^2./N(idm).^2)));
%dV(idp) = 0;
%dV(idm) = 0;
