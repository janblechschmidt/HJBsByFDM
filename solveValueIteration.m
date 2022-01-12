% TODO: The following is an excerpt and should be completed
% x_prev = w_prev;
% VIdone = 0;
% while ~VIdone
% 	dx = S*x_prev - rhs;
% 	x = x_prev - 0.002 * dx;
% 	VI_norm = a.norm1(x - x_prev);
% 	% sfigure(1); clf;
% 	% plotValue(p, c, a, x)
% 	% pause
% 	if VI_norm < 1e-8
% 		VIdone = 1;
% 	else
% 		x_prev = x;
% 	end
% end
% w = x;
