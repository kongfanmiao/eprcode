function resM = unitary_transform(M,U)
% Do unitary transformation
% M: orginal matrix
% U: unitary matrix

% check if U is unitary
Unorm = U*U';
if isa(Unorm,'sym')
    Unorm = simplify(Unorm);
end
Unorm = abs(double(Unorm));
if ~all(ismembertol(Unorm,eye(size(U))),'all')
    error("Please provide a unitary matrix");
end
resM = U*M*U';
end