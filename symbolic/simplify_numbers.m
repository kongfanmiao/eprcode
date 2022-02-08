function res = simplify_numbers(x)
% remove tiny numbers, tolerance is 1e-8
if abs(x) <= 1e-8
    tmp = 0;
else
    tmp = x;
end 
% simplify the fractional numbers
[N,D] = rat(tmp);
res = N./D;
end