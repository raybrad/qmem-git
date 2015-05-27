function scaling_static()

global epsilon q scl;
global T;
global nodes nodeV linkL linkS dlinkL linkVolumes nodeVolumes;

%%% three independent scaling parameters (lambda is changeable),default 1e-9
scl = struct('T',300,'Vt',2.5852,'lambda',1e-9,'ni',[],'s_D',[],'s_mu',[],...
    's_J',[],'tao',[],'s_E',[],'omega',[],'s_sigma',[],'s_Curr',[],...
    's_v',[],'s_A',[],'K',[]);

scl.ni = epsilon*scl.Vt/q/scl.lambda^2; % C/V/m * V / C / m^2  ~ 1/m^3 
scl.s_E = scl.Vt/scl.lambda;
T = T/scl.T;
%the 1e-3 is deleted from initGeometry,here 1e-9 is used.
nodes = nodes*1e-9/scl.lambda;
nodeV = nodeV*1e-27./(scl.lambda^3);
linkL = linkL*1e-9./scl.lambda;
linkS = linkS*1e-18./(scl.lambda^2);
dlinkL = dlinkL*1e-9./scl.lambda;
for i = 1:length(linkVolumes)
   linkVolumes{i}(2,:) = linkVolumes{i}(2,:)*1e-18/(scl.lambda^2);
end
for i = 1:length(nodeVolumes)
   nodeVolumes{i}(2,:) = nodeVolumes{i}(2,:)*1e-27/(scl.lambda^3);
end

