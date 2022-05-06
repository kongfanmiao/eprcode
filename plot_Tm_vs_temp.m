function plot_Tm_vs_temp(FitResult,varargin)

par = inputParser;
par.addOptional('inverseT1', '', @(x)any(validatestring(x, ...
    {'', 'inverse'})));
par.addParameter('Color', 'b', @(x)isstring(x)|ischar(x));
par.addParameter('Marker', 'o', @(x)any(validatestring(x, ...
    {'o','+','*','.','x','-','|','s','d','^','v','>','<','p','h'})));
par.addParameter('MarkerSize', 100, @isnumeric);
par.KeepUnmatched = true;

parse(par, varargin{:});
args = par.Results;

temp = FitResult.Temperature;

try
    Tm = FitResult.Tm; % mono or stretch exponential
catch
    try
        Tm = FitResult.T_long; % bi-exponential
    catch
        error('There is no Tm or T_long in the FitResult table');
    end
end

if ~isempty(args.inverseT1)
    y = 1./Tm;
    ylab = sprintf('{T_m}^{-1} (%s^{-1})', ...
        FitResult.Properties.VariableUnits{2});
    ttl = '{T_m}^{-1} versus Temperature';
else
    y = Tm;
    ylab = sprintf('T_m (%s)', ...
        FitResult.Properties.VariableUnits{2});
    ttl = 'T_m versus Temperature';
end

set(gcf,'color','w');
scatter(temp, y, args.MarkerSize, args.Marker, args.Color,'filled')
xlabel("Temperature (K)");
ylabel(ylab);
title(ttl)
box on

fig = gcf;
if fig.WindowStyle ~= "docked"
    set(fig,'position',[10,10,900,600]);
end


end