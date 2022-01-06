function rounded_dm = roundDM(dm,tol)
% round the double elements in density matrix
% default tolerance is 8
arguments
    dm
    tol = 8
end
for i = 1:numel(dm)
    el = dm(i);
    try
        dm(i) = round(double(el),tol);
    catch
    end
end
rounded_dm = dm;
end