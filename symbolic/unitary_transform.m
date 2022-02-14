function resM = unitary_transform(M,U)
% Do unitary transformation
% M: orginal matrix
% U: unitary matrix

% % check if U is unitary
% Unorm = U*U';
% if isa(Unorm,'sym')
%     Unorm = simplify(Unorm);
%     Unorm = expand(Unorm,"ArithmeticOnly", true);
%     Unorm = simplify(Unorm);
% end
% Unorm = abs(double(Unorm));
% if ~all(ismembertol(Unorm,eye(size(U))),'all')
%     error("Please provide a unitary matrix");
% end
resM = U*M*U';
if isa(resM, 'sym')
    resM = mapSymType(resM, 'rational', @(x)simplify_numbers(x));
    resM = expand(resM, "ArithmeticOnly", true);
    resM = combine(resM, "sincos");
end
end