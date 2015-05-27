function [boundary,internal] = defBoundaryNodes(numx,numy,numz)

%**********************************************************************%
%----------------------------------------------------------------------%
%|    This function defines the grid boundaryry                       |%
%|                                                                    |%
%|                                                                    |%  
%|    Civil and Computational Engineering Centre                      |%       
%|    University of Wales Swansea                                     |%  
%|                                                                    |%  
%|                                                                    |%  
%----------------------------------------------------------------------%
%**********************************************************************%

%
% Synopsis:   [boundary,internal,domain] = def_boundary(numx,numy,numz)

%
% Input:  
%         numx              >  Number of control volumes in X direction.
%         numy              >  Number of control volumes in Y direction.
%         numz              >  Number of control volumes in Z direction.
%
% Output: 
%         boundary          >  Boundary Elements.
%         internal          >  Internal Elements of the grid.
%         domain            >  Domain of the grid


% The boundary elements of the entire domain are established %%%%%%

domain = 1:numx*numy*numz;
bottom_boundary  = 1:numx*numy;
top_boundary    = (numx*numy*(numz-1))+1:numx*numy*numz;
left_boundary   = 1:numx:(numx*numy*numz-numx+1);
right_boundary  = numx:numx:numx*numy*numz;

for k = 1:numy
    for j = 1:numz
        for i = 1:numx
            dom(i,j,k) = i+ numx*(k-1) + (j-1)*numx*numy;
        end
    end
end
front_boundary    = reshape(dom(:,:,1),numx*numz,1)'; 
back_boundary = reshape(dom(:,:,numy),numx*numz,1)'; 


boundary=[front_boundary back_boundary left_boundary right_boundary bottom_boundary top_boundary ]; %Boundary is comprised of all the edge elements% 
boundary= unique(boundary);
  
internal = setxor(domain,boundary);% Internal cells are all the cells that 
                                  % are not part of the boundary
boundary = boundary';
internal = internal';
domain = domain';

