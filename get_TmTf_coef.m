function result = get_TmTf_coef(FitResult)
% Extract Tm and Tf and their coefficients from bi-exponential fitting
arguments
    FitResult table
end

temperature = FitResult.temperature;
TmTf = 1./FitResult.k;
[Tm, TmIdx] = max(TmTf,[],2);
[Tf, TfIdx] = min(TmTf,[],2);
TmCoef = zeros(size(TmIdx));
TfCoef = zeros(size(TfIdx));
for i = 1:length(TmIdx)
    TmCoef(i) = -FitResult.c(i,TmIdx(i)+1);
    TfCoef(i) = -FitResult.c(i,TfIdx(i)+1);
end
result = table(temperature, Tm, TmCoef, Tf, TfCoef);
end