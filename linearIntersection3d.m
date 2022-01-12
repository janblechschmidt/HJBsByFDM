function Xnew = linearIntersection3d(X, alpha, c)
% This function is intended to compute the points in a tetrahedron,
% where a condition alpha(x) = c is attained with equality.
% Typically, there are alpha(x) > c or alpha(x) < c.
% For tetrahedrons where both cases occur, a boundary alpha(x) = c should lie
% within the boundaries of this tetrahedron.
% If we assume a linear connection, the intersection is either a triangle or
% a quadrangle (if the case is not degenerated).

abool = (alpha <= c);
[n,m] = size(X);
Xnew = zeros(0,m);

if any(abool) ~= all(abool)
	k = 1;
	switch mod(sum(abool),2)
	case 0
		% Quadrilateral
		% Sorting is important, otherwise matlab will plot two triangles
		if allOrNothing(abool', [1, 0, 0, 1])
			i = 1; j = 2;
			l = (c - alpha(i)) / (alpha(j) - alpha(i));
			Xnew(1,:) = (1.-l) * X(i,:) + l * X(j,:);
			i = 2; j = 4;
			l = (c - alpha(i)) / (alpha(j) - alpha(i));
			Xnew(2,:) = (1.-l) * X(i,:) + l * X(j,:);
			i = 4; j = 3;
			l = (c - alpha(i)) / (alpha(j) - alpha(i));
			Xnew(3,:) = (1.-l) * X(i,:) + l * X(j,:);
			i = 3; j = 1;
			l = (c - alpha(i)) / (alpha(j) - alpha(i));
			Xnew(4,:) = (1.-l) * X(i,:) + l * X(j,:);
		else
			if allOrNothing(abool', [1, 1, 0, 0])
				i = 1; j = 3;
				l = (c - alpha(i)) / (alpha(j) - alpha(i));
				Xnew(1,:) = (1.-l) * X(i,:) + l * X(j,:);
				i = 3; j = 2;
				l = (c - alpha(i)) / (alpha(j) - alpha(i));
				Xnew(2,:) = (1.-l) * X(i,:) + l * X(j,:);
				i = 2; j = 4;
				l = (c - alpha(i)) / (alpha(j) - alpha(i));
				Xnew(3,:) = (1.-l) * X(i,:) + l * X(j,:);
				i = 4; j = 1;
				l = (c - alpha(i)) / (alpha(j) - alpha(i));
				Xnew(4,:) = (1.-l) * X(i,:) + l * X(j,:);
			else
				if allOrNothing(abool', [1, 0, 1, 0])
					i = 1; j = 2;
					l = (c - alpha(i)) / (alpha(j) - alpha(i));
					Xnew(1,:) = (1.-l) * X(i,:) + l * X(j,:);
					i = 2; j = 3;
					l = (c - alpha(i)) / (alpha(j) - alpha(i));
					Xnew(2,:) = (1.-l) * X(i,:) + l * X(j,:);
					i = 3; j = 4;
					l = (c - alpha(i)) / (alpha(j) - alpha(i));
					Xnew(3,:) = (1.-l) * X(i,:) + l * X(j,:);
					i = 4; j = 1;
					l = (c - alpha(i)) / (alpha(j) - alpha(i));
					Xnew(4,:) = (1.-l) * X(i,:) + l * X(j,:);
				else
					error('no case\n')
				end % if allOrNothing(abool', [1, 0, 1, 0])
			end % if allOrNothing(abool', [1, 1, 0, 0])
		end %if allOrNothing(abool', [1, 0, 0, 1])
	case 1
		% Triangle
		for i = 1:m
			for j = i+1:m+1
				if abool(i) ~= abool(j)
					l = (c - alpha(i)) / (alpha(j) - alpha(i));
					Xnew(k,:) = (1.-l) * X(i,:) + l * X(j,:);
					k = k + 1;
				end % if abool(i) ~= abool(j)
			end % for j = i+1:m+1
		end % for i = 1:m
	end % switch mod(sum(abool),2)
end % if any(abool) != all(abool):

	% 	if np.mod(sum(abool),2) == 1:
	% se:
	% 	# Intersection is quadrangle
	% 	Xnew = np.ndarray((dim+1,dim),dtype=np.double)
	% 	if allOrNothing(abool == np.array([True, False, False, True])):
	% 			# print('case 1')
	% 			i = 0; j = 1
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[0] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			i = 1; j = 3
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[1] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			i = 3; j = 2
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[2] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			i = 2; j = 0
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[3] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			elif allOrNothing(abool == np.array([True, True, False, False])):
	% 			# print('case 2')
	% 			i = 0; j = 2
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[0] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			i = 2; j = 1
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[1] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			i = 1; j = 3
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[2] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			i = 3; j = 0
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[3] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			elif allOrNothing(abool == np.array([True, False, True, False])):
	% 			# print('case 3')
	% 			i = 0; j = 1
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[0] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			i = 1; j = 2
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[1] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			i = 2; j = 3
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[2] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 			i = 3; j = 0
	% 			l = (c - alpha[i]) / (alpha[j] - alpha[i])
	% 			Xnew[3] = (1.-l) * Xalt[i] + l * Xalt[j]
	% 	else:
	% 			print('unknown case')
	% 			sys.exit()
	% 			return Xnew
	% 	else:
	% 			return False
