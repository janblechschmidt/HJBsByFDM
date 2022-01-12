function A = setupDifferentialOperator(p, a, mode, dim, dim2, mixed_mode)
A = sparse(p.ndofs,p.ndofs);
if nargin < 5
	dim2 = -1;
end

if strcmp(mode, 'FirstOrder')
	if p.dim == 1
		switch dim
		case 1
			idx_central  = setdiff(a.iall, unique([a.ixlow; a.ixup]));
			idx_forward  = [a.ixlow];
			idx_backward = [a.ixup];

        end % switch dim

    end % if p.dim == 1

	if p.dim == 2
		switch dim
		case 1
			% idx_central  = [a.iinner; a.iylowinner; a.iyupinner];
			idx_central  = setdiff(a.iall, unique([a.ixlow; a.ixup]));
			idx_forward  = [a.ixlow];
			idx_backward = [a.ixup];

		case 2
			% idx_central  = [a.iinner; a.ixlowinner; a.ixupinner];
			idx_central  = setdiff(a.iall, unique([a.iylow; a.iyup]));
			idx_forward  = [a.iylow];
			idx_backward = [a.iyup];
		end % switch dim
	end % if p.dim == 2

	if p.dim == 3
		switch dim
		case 1
			idx_central  = setdiff(a.iall, unique([a.ixlow; a.ixup]));
			idx_forward  = [a.ixlow];
			idx_backward = [a.ixup];

		case 2
			idx_central  = setdiff(a.iall, unique([a.iylow; a.iyup]));
			idx_forward  = [a.iylow];
			idx_backward = [a.iyup];

		case 3
			idx_central  = setdiff(a.iall, unique([a.izlow; a.izup]));
			idx_forward  = [a.izlow];
			idx_backward = [a.izup];

		end % switch dim
	end % if p.dim == 3
end % if strcmp(mode='FirstOrder')

if strcmp(mode,'FirstOrderForward')
	if p.dim == 1
		switch dim
		case 1
			idx_central  = [];
			idx_forward  = setdiff(a.iall, a.ixup);
			idx_backward = [a.ixup];

		end % switch dim
	end % if p.dim == 1

	if p.dim == 2
		switch dim
		case 1
			idx_central  = [];
			idx_forward  = setdiff(a.iall, a.ixup);
			% idx_forward  = [a.iinner; a.ixlow; a.iylowinner; a.iyupinner];
			idx_backward = [a.ixup];

		case 2
			idx_central  = [];
			idx_forward  = setdiff(a.iall, a.iyup);
			% idx_forward  = [a.iinner; a.iylow; a.ixlowinner; a.ixupinner];
			idx_backward = [a.iyup];
		end % switch dim
	end % if p.dim == 2

	if p.dim == 3
		switch dim
		case 1
			idx_central  = [];
			idx_forward  = setdiff(a.iall, a.ixup);
			% idx_forward  = unique([a.iinner; a.ixlow; a.iylowinner; a.iyupinner; a.izlowinner; a.izupinner]);
			idx_backward = [a.ixup];

		case 2
			idx_central  = [];
			idx_forward  = setdiff(a.iall, a.iyup);
			% idx_forward  = unique([a.iinner; a.iylow; a.ixlowinner; a.ixupinner; a.izlowinner; a.izupinner]);
			idx_backward = [a.iyup];

		case 3
			idx_central  = [];
			idx_forward  = setdiff(a.iall, a.izup);
			% idx_forward  = unique([a.iinner; a.izlow; a.ixlowinner; a.ixupinner; a.iylowinner; a.iyupinner]);
			idx_backward = [a.izup];
		end % switch dim

	end % if p.dim == 3
end % if strcmp(mode,'FirstOrderForward')

if strcmp(mode,'FirstOrderBackward')
	if p.dim == 1
		switch dim
		case 1
			idx_central  = [];
			idx_forward  = [a.ixlow];
			idx_backward  = setdiff(a.iall, a.ixlow);
		end % switch dim
	end % if p.dim == 1

	if p.dim == 2
		switch dim
		case 1
			idx_central  = [];
			idx_forward  = [a.ixlow];
			idx_backward  = setdiff(a.iall, a.ixlow);
			% idx_backward  = [a.iinner; a.ixup; a.iylowinner; a.iyupinner];

		case 2
			idx_central  = [];
			idx_forward  = [a.iylow];
			idx_backward  = setdiff(a.iall, a.iylow);
			% idx_backward  = [a.iinner; a.iyup; a.ixlowinner; a.ixupinner];
		end % switch dim
	end % if p.dim == 2

	if p.dim == 3
		switch dim
		case 1
			idx_central  = [];
			idx_forward  = [a.ixlow];
			idx_backward  = setdiff(a.iall, a.ixlow);
			% idx_backward  = unique([a.iinner; a.ixup; a.iylowinner; a.iyupinner; a.izlowinner; a.izupinner]);

		case 2
			idx_central  = [];
			idx_forward  = [a.iylow];
			idx_backward  = setdiff(a.iall, a.iylow);
			% idx_backward  = unique([a.iinner; a.iyup; a.ixlowinner; a.ixupinner; a.izlowinner; a.izupinner]);

		case 3
			idx_central  = [];
			idx_forward  = [a.izlow];
			idx_backward  = setdiff(a.iall, a.izlow);
			% idx_backward  = unique([a.iinner; a.izup; a.ixlowinner; a.ixupinner; a.iylowinner; a.iyupinner]);

		end % switch dim
	end % if p.dim == 3
end % if strcmp(mode,'FirstOrderBackward')

% The following has to be done for all first order approximations
if strfind(mode, 'First')
	% TODO: Sicherung, dass summe(1) == size(a.XY,0) ist
	if size(idx_central,1) + size(idx_forward,1) + size(idx_backward,1) ~= size(a.XY,1)
		warning('not all points captured by first order differential operator')
		keyboard
	end
	
	% Central difference for idx_central
	if ~isempty(idx_central)
			[idx, vals] = stencil_FirstOrder(a, 'Central', idx_central, dim);
			for i=1:size(idx,2)
				A = A + sparse(idx_central, idx(:,i), vals(:,i), p.ndofs, p.ndofs);
			end
		end % if ~isempty(idx_central)
	
		% Central difference for idx_forward
		if ~isempty(idx_forward)
			[idx, vals] = stencil_FirstOrder(a, 'Forward', idx_forward, dim);
			for i=1:size(idx,2)
				A = A + sparse(idx_forward, idx(:,i), vals(:,i), p.ndofs, p.ndofs);
			end
		end % if ~isempty(idx_forward)
	
		% Central difference for idx_backward
		if ~isempty(idx_backward)
			[idx, vals] = stencil_FirstOrder(a, 'Backward', idx_backward, dim);
			for i=1:size(idx,2)
				A = A + sparse(idx_backward, idx(:,i), vals(:,i), p.ndofs, p.ndofs);
			end
		end % if ~isempty(idx_backward)
end % if strfind(mode, 'First')

if strcmp(mode,'SecondOrder')
	if p.dim == 1
		switch dim
		case 1
			idx_second = a.ixinner;
		end % switch dim
	end % if p.dim == 1

	if p.dim == 2
		switch dim
		case 1
			idx_second = a.ixinner;
		case 2
			idx_second = a.iyinner;
		end % switch dim
	end % if p.dim == 2

	if p.dim == 3
		switch dim
		case 1
			idx_second = a.ixinner;
		case 2
			idx_second = a.iyinner;
		case 3
			idx_second = a.izinner;
		end % switch dim
	end % if p.dim == 3

	[idx, vals] = stencil_SecondOrder(a, idx_second, dim);
	for i=1:size(idx,2)
		A = A + sparse(idx_second, idx(:,i), vals(:,i), p.ndofs, p.ndofs);
	end
	% keyboard
end % if strcmp(mode,'SecondOrder')

if strcmp(mode,'Mixed')

	if p.dim == 2
		% Only x and y possible
		idx_second = [a.iinner];
	end % if p.dim == 2

	if p.dim == 3
		switch setdiff(1:3,[dim, dim2])
		case 1
			idx_second = [a.iinner; a.ixlowinner; a.ixupinner];
		case 2
			idx_second = [a.iinner; a.iylowinner; a.iyupinner];
		case 3
			idx_second = [a.iinner; a.izlowinner; a.izupinner];
		end % switch dim
	end % if p.dim == 3

	[idx, vals] = stencil_SecondMixed(a, idx_second, mixed_mode, dim, dim2);
	for i=1:size(idx,2)
		A = A + sparse(idx_second, idx(:,i), vals(:,i), p.ndofs, p.ndofs);
	end
end % if strcmp(mode,'SecondOrder')

