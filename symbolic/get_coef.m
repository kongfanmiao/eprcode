function coef = get_coef(op,base,digits,tol)
arguments
    op
    base
    digits = 8 
    tol = 1e-8
end
% constant = sqrt(trace(base*base'));
% if constant ~= 1
%     error("The base should be normalized")
% end
coef = trace(base*op);
if isa(coef, 'sym')
    coef = simplify(coef);
    coef = combine(coef, 'sincos');
    % remove tiny numbers
    coef = mapSymType(coef, 'rational', ...
        @(x)piecewise(abs(x)<=tol, 0, x));
    % express fraction in floating point
%     coef = sym(vpa(coef, digits));
elseif isa(coef, 'float')
    coef = round(coef, digits);
end
end
