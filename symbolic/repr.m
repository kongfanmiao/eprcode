function [op_repr, op_items] = repr(op)

load("operators.mat");
% orthonormal basis
% basis_str = {'I','S_x','S_y','S_z','I_x','I_y','I_z', ...
%     'S_z*I_z','S_x*I_z','S_y*I_z','S_z*I_x','S_z*I_y', ...
%     'S_x*I_x','S_y*I_y','S_x*I_y','S_y*I_x'};
basis_str = {'I','Sx','Sy','Sz','Ix','Iy','Iz', ...
    'SzIz','SxIz','SyIz','SzIx','SzIy', ...
    'SxIx','SyIy','SxIy','SyIx'};
op_repr = 0;
op_items.basis = {};
op_items.coeffients = {};
% express in linear combinatino of basis
for i = 1:16
    b = basis{i};
    constant = sqrt(trace(b*b'));
    coef = get_coef(op,b/constant);
    b_str = basis_str{i};
    b_sym = sym(b_str);
    if coef ~= 0
        op_items.basis{end+1} = b;
        op_items.coeffients{end+1} = coef;
    end
    op_repr = op_repr + coef*b_sym/constant;
end
% simplify the expression
if isa(op_repr, 'sym')
    op_repr = simplify(op_repr);
end
end