function [Tl, TlPerc, Ts, TsPerc] = ...
    get_Tlong_and_Tshort(k, c, timeScaleFactor)
% Extract T_long, T_short and their coefficients from bi-exponential fitting

arguments
    k (:,2) double
    c (:,3) double
    timeScaleFactor double
end

TlTs = k/timeScaleFactor;
[Tl, TlIdx] = max(TlTs,[],2); 
[Ts, TsIdx] = min(TlTs,[],2);
TlCoef = zeros(size(TlIdx));
TsCoef = zeros(size(TsIdx));
for i = 1:length(TlIdx)
    TlCoef(i) = c(i,TlIdx(i)+1);
    TsCoef(i) = c(i,TsIdx(i)+1);
end
TlPerc = abs(TlCoef)./(abs(TlCoef)+abs(TsCoef));
TsPerc = abs(TsCoef)./(abs(TlCoef)+abs(TsCoef));
end