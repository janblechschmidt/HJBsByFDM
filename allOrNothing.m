function res = allOrNothing(a,b)
% allOrNothing compares to arrays of logical expressions and yiels true,
% if all elements are the same or all elements are differing.
if all(a == b) || all(a ~= b)
	res = 1;
else 
	res = 0;
end
