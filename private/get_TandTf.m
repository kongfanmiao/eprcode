function [T, TPerc, Tf, TfPerc] = get_TandTf(k, c, timeScaleFactor)
% Extract T and Tf and their coefficients from bi-exponential fitting

arguments
    k (:,2) double
    c (:,3) double
    timeScaleFactor double
end

TTf = k/timeScaleFactor;
% The larger contribution is made by T, smaller by Tf
[T, TIdx] = max(TTf,[],2);
[Tf, TfIdx] = min(TTf,[],2);
TCoef = zeros(size(TIdx));
TfCoef = zeros(size(TfIdx));
for i = 1:length(TIdx)
    TCoef(i) = c(i,TIdx(i)+1);
    TfCoef(i) = c(i,TfIdx(i)+1);
end
TPerc = abs(TCoef)./(abs(TCoef)+abs(TfCoef));
TfPerc = abs(TfCoef)./(abs(TCoef)+abs(TfCoef));
end