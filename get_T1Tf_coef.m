function result = get_T1Tf_coef(FitResult)
% Extract T1 and Tf and their coefficients from bi-exponential fitting
arguments
    FitResult table
end

temperature = FitResult.temperature;
T1Tf = 1./FitResult.k;
[T1, T1Idx] = max(T1Tf,[],2);
[Tf, TfIdx] = min(T1Tf,[],2);
T1Coef = zeros(size(T1Idx));
TfCoef = zeros(size(TfIdx));
for i = 1:length(T1Idx)
    T1Coef(i) = -FitResult.c(i,T1Idx(i)+1);
    TfCoef(i) = -FitResult.c(i,TfIdx(i)+1);
end
result = table(temperature, T1, T1Coef, Tf, TfCoef);
end
    

