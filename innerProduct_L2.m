function ip = innerProduct_L2(h, v, w)
% innerProduct_L2 computes the discrete L^2 inner product of two functions v, w
% INPUT: h - array of equidistant distance between two points in spatial dimensions
%        v - first discrete function on function space V_h
%        w - second discrete function on function space V_h

ip = prod(h) * sum(v.*w);

end % function ip = innerProduct_L2(h, v, w)
