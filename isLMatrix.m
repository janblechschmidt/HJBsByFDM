function val = isLMatrix(A)
% Function to prove whether matrix A is a L-matrix
% i.e. all off-diagonal entries are non-positive and the diagonal entries are positive
% i.e. A_i,j <= 0 for all i=j and A_i,i > 0 for all i

AdiagZero = A - diag(diag(A));

if max(max(AdiagZero)) <= 0 & min(min(diag(A))) > 0
	val = 1;
else
	val = 0;
end

end % function isZMatrix(A)
