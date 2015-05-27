function scaling_dynamic(omega1)

global epsilon sigma0 sigma mu q scl;
global omega v omega_p gamma_p;
global lomega_1 lgamma_1 lomega_2 lgamma_2; %Lorentz pole
global dt epdf epdf2 light_speed;
global tlas tzero;

scl.s_D = 1;	%  m^2/s
scl.tao = scl.lambda^2/scl.s_D;

scl.s_mu = scl.s_D/scl.Vt;
scl.omega = 1/scl.tao;
scl.s_sigma = epsilon/scl.tao;
scl.s_A = scl.tao*scl.Vt/scl.lambda;
scl.s_B = scl.tao*scl.Vt/scl.lambda^2;
scl.s_E = scl.Vt/scl.lambda;
scl.s_J = q*scl.ni*scl.s_D/scl.lambda;	% C /m^3 * m^2/s /m ~  C/s/m^2 ~ A/m^2  %epsilon*Vt/lambda^3 or epsilon*Vt/lambda/tao
scl.s_Curr = scl.s_sigma * scl.Vt * scl.lambda;
scl.K = epsilon*mu*(scl.lambda/scl.tao)^2;
v = sigma/(omega*epsilon);
sigma = sigma0/scl.s_sigma;
omega = omega/scl.omega;
dt  = dt/scl.tao;
tlas  = tlas/scl.tao;
tzero  = tzero/scl.tao;
epdf = epdf/scl.tao;

%light source & plasmonic metal
epdf2 = epdf2/scl.tao;
light_speed=sqrt(1/(epsilon*mu))/(scl.lambda/scl.tao);
omega_p=omega_p/scl.omega;
gamma_p=gamma_p/scl.omega;

lomega_1=lomega_1/scl.omega;
lgamma_1=lgamma_1/scl.omega;
lomega_2=lomega_2/scl.omega;
lgamma_2=lgamma_2/scl.omega;
