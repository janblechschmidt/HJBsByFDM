function val = isStrictlyDiagonallyDominant(A)
% Function to prove whether matrix A is a strictly diagonally dominant which implies the M-matrix property

diagvals = abs(diag(A));
AdiagZero = A - diag(diag(A));
offdiagvals = sum(abs(AdiagZero),2);

if all(diagvals > offdiagvals)
	val = 1;
else
	val = 0;
end % if all(diagvals > offdiagvals)
