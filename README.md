Version 1.0
%Revised version 1
%0 light source is added as a sheet of current density Js, which create plane wave field 
%0 Mur absorbing boundary is added to reduce reflection
%2 plasmonic metal is modelled using Two Lorentz poles and one drude pole.
%3 tdrelaxstep2.m use (row,col,val) technique to create sparse matrix
%3 in initGeometry sparse matrix creation is used for linkVolS nodeVolV

%Revised version 2
%0 light source is added as a sheet of current density Js, which create plane wave field 
%0 Mur absorbing boundary is added to reduce reflection
%2 plasmonic metal is modelled using Two Lorentz poles and one drude pole.
%3 tdrelaxstep2.m use (row,col,val) technique to create sparse matrix
%3 in initGeometry sparse matrix creation is used for linkVolS nodeVolV
%4 tdbuilJacob.m prebuild,precondition and save  the Jacobian matrix,which is unchanged in time propagation. Semiconductor part is removed.
%4 tdsolverhs construct the right hand side matrix, and some simple vectorization technique is used. So tdrelaxstep2 is replaced.

%Revised version 3
%0 light source is added as a sheet of current density Js, which create plane wave field 
%0 Mur absorbing boundary is added to reduce reflection
%2 plasmonic metal is modelled using Two Lorentz poles and one drude pole.
%3 tdrelaxstep2.m use (row,col,val) technique to create sparse matrix
%3 in initGeometry sparse matrix creation is used for linkVolS nodeVolV
%4 tdbuilJacob.m prebuild,precondition and save  the Jacobian matrix,which is unchanged in time propagation. Semiconductor part is removed.
%4 tdsolverhs construct the right hand side matrix, and some simple vectorization technique is used. So tdrelaxstep2 is replaced.
%5 tdbuildRHSCoef fully vectorize the rhs coefficients,prebuild and save it. 
%5 tdsolerhs_vector use previous coefficient matrix to reconstruct the rhs matrix. So tdsolverhs is replaced. 
%6 in initGeometry the allocation of dirlinks is added to speed up
%6 in scaling_static change Vt from 2.5852e-2 to 2.5852,so to rescale a little bit to balance the num of Js by scl.J
%6 in initGeometry linkVolumes nodeVolumes are deleted and just used linkVolS and nodeVolV are ok,variables are changed in following sub
