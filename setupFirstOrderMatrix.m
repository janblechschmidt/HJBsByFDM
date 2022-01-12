function A1func  = setupFirstOrderMatrix(p, a)
% setupFirstOrderMatrix sets up the matrix for all first-order derivative terms
% defined in a problem structure p and an algorithmic structure a.
% This is done using the convection terms in the coefficient strucure

% Load helper functions
helperFunctions;

% Initialize dimension of first order matrix
m = p.ndofs;
n = p.ndofs;

% Evaluate convection coeffcients at grid points
% for i = 1:p.dim
%     conv{i} = @(t,u) p.convection{i}(a.XY, t, u);
% end

if ~isfield(p, 'deterministicDimension')
    p.deterministicDimension = [];
end

A1func = @(t,u) sparse(m,n);
if p.useUpwinding

    % New implementation

	for i = 1:p.dim
        if ~ismember(i, p.deterministicDimension) || strcmp(p.solver, 'Howard')
            % fprintf('First order setup upwind dim %d\n', i);
            Aforw = setupDifferentialOperator(p, a, 'FirstOrderForward', i);
    		Aback = setupDifferentialOperator(p, a, 'FirstOrderBackward', i);

            % Evaluate convection coeffcients at grid points
            conv = @(t,u) p.convection{i}(a.XY, t, u);

            A1func = @(t, u) A1func(t, u) + ...
                spdiags(posPart(conv(t,u)),0,m,n) * Aforw - ...
                spdiags(negPart(conv(t,u)),0,m,n) * Aback;
        end
	end % for i = 1:p.dim


	% % Get matrices for differential operators
	% for i = 1:p.dim
	%     Aforw{i} = setupDifferentialOperator(p, a, 'FirstOrderForward', i);
	%     Aback{i} = setupDifferentialOperator(p, a, 'FirstOrderBackward', i);
	% end % for i = 1:p.dim

	% switch p.dim
	% case 1
	%     A1func = @(t, u) spdiags(posPart(conv{1}(t,u)),0,m,n) * Aforw{1}  - ...
	%             spdiags(negPart(conv{1}(t,u)),0,m,n) * Aback{1};

	% case 2
    %     A1func = @(t, u) spdiags(posPart(conv{1}(t,u)),0,m,n) * Aforw{1}  - ...
    %             spdiags(negPart(conv{1}(t,u)),0,m,n) * Aback{1} + ...
    %             spdiags(posPart(conv{2}(t,u)),0,m,n) * Aforw{2} - ...
    %             spdiags(negPart(conv{2}(t,u)),0,m,n) * Aback{2};

	% case 3
	%     A1func = @(t, u) spdiags(posPart(conv{1}(t,u)),0,m,n) * Aforw{1}  - ...
	%             spdiags(negPart(conv{1}(t,u)),0,m,n) * Aback{1} + ...
	%             spdiags(posPart(conv{2}(t,u)),0,m,n) * Aforw{2} - ...
	%             spdiags(negPart(conv{2}(t,u)),0,m,n) * Aback{2} + ...
	%             spdiags(posPart(conv{3}(t,u)),0,m,n) * Aforw{3} - ...
	%             spdiags(negPart(conv{3}(t,u)),0,m,n) * Aback{3};
	% end % switch p.dim

else

    % New implementation
	for i = 1:p.dim
        if ~ismember(i,p.deterministicDimension) || strcmp(p.solver, 'Howard')
            % fprintf('First order setup central dim %d', i);
		    Ac = setupDifferentialOperator(p, a, 'FirstOrder', i);

            % Evaluate convection coeffcients at grid points
            conv = @(t,u) p.convection{i}(a.XY, t, u);

		    A1func = @(t, u) A1func(t, u) + spdiags(conv(t,u), 0, m, n) * Ac;
        end
	end % for i = 1:p.dim

	% % Get matrices for differential operators
	% for i = 1:p.dim
	%     A{i} = setupDifferentialOperator(p, a, 'FirstOrder', i);
	% end % for i = 1:p.dim

	% switch p.dim
	% case 1
	%     A1func = @(t, u) spdiags(conv{1}(t,u),0,m,n) * A{1};

	% case 2
	%     A1func = @(t, u) spdiags(conv{1}(t,u),0,m,n) * A{1} + ...
	%             spdiags(conv{2}(t,u),0,m,n) * A{2};

	% case 3
	%     A1func = @(t, u) spdiags(conv{1}(t,u),0,m,n) * A1{1} + ...
	%             spdiags(conv{2}(t,u),0,m,n) * A{2};
	%             spdiags(conv{3}(t,u),0,m,n) * A{3};
	% end % switch p.dim

end % if p.useUpwinding

end % function A1func  = setupFirstOrderMatrix(p, a)
