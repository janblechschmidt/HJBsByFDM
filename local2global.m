% This function converts the index vectors ix1, ix2 (ix3) referring to 
% x1, x2 (and x3) node numbers (ranging from 0 to Nx1, or 0 to Nx2, (or 0
% to Nx3) respectively), into one vector of global indices 
% (ranging from 1 to (Nx1+1)*(Nx2+1)*(Nx3+1)), according to 
% the lexicographic ordering.
% ix1 is the most significant index.
% The program can run with 4 (2d version) or 6 (3d version) input parameters.

function ix = local2global(ix1,Nx1,ix2,Nx2,ix3,Nx3)

if nargin == 4
	ix = (ix2 + 1) + ix1 * (Nx2+1);
else
	if nargin == 6
		ix = (ix3 + 1) + ix2 * (Nx3+1) + ix1 * (Nx2+1)*(Nx3+1);
	else
		error('Unknown local-to-global mapping.')
	end % if nargin == 6
end % if nargin == 4

end % function ix = local2global(ix1,Nx1,ix2,Nx2,ix3,Nx3)
