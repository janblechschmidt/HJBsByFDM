function plotDecisionBounday3d(decBnd, type, varargin)

if strcmpi(type, 'triangle')
	d = 3;
end % if strcmpi(type, 'triangle')

if strcmpi(type, 'quadrilateral')
	d = 4;
end % if strcmpi(type, 'quadrilateral')

n = size(decBnd,1)/d;

X = reshape(decBnd(:,1),d,n);
Y = reshape(decBnd(:,2),d,n);
Z = reshape(decBnd(:,3),d,n);
% fill3(X, Y, Z, varargin{:})
patch(X, Y, Z, varargin{:})
axis([0,1,0,1,0,1])
