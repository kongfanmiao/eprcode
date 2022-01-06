function result = spin_evolve(sigma, Hamiltonian)
% express B in Cartesian basis
% sigma: initial density matrix
% Hamiltonian: Hamiltonian that operats on the density matrix, should
% include the time

load('operators.mat');

result = sigma;
prev = 1;
% check if the terms in Hamiltonian commute with each other
for i = 1:16
    b = basis{i};
    constant = sqrt(trace(b*b'));
    coef = get_coef(Hamiltonian,b/constant);
    if coef ~= 0
        if ~all(commutator(b,prev)==0)
            error("The terms in Hamiltonian should commute with each other");
        end
%         result = result*cos(coef) + ...
%             1i*commutator(result,b/constant)*sin(coef);
        prev = b;
    end
end
U = expm(-1i*Hamiltonian);
result = unitary_transform(sigma,U);
end