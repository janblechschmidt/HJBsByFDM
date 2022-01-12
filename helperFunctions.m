% This script defines some nice helper functions

% posPart gives the positive part of a scalar, vector or matrix in the same format
posPart = @(x) x .* (x > 0);

% negPart gives the negative part of a scalar, vector or matrix in the same format
negPart = @(x) -x .* (x < 0);

% zeroDiag returns the matrix with zero diagonal entries
zeroDiag = @(A) A - diag(diag(A));
