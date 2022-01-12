function n = norm_L2(h, v)
% Computes the discrete L^2 norm of a function v
% INPUT: h - array of equidistant distance between two points in spatial dimensions
%        v - discrete function on function space V_h
% n = prod(h) * sum(v.^2);
n = sum(h.*(v.^2));

end % function n = norm_L2(h, v)
