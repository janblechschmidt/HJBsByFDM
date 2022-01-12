function val = isZMatrix(A)
% Function to prove whether matrix A is a Z-matrix
% i.e. all off-diagonal entries are non-positive
% i.e. A_i,j <= 0 for all i=j

AdiagZero = A - diag(diag(A));

if max(max(AdiagZero)) <= 0
	val = 1;
else
	val = 0;
end

end % function isZMatrix(A)
