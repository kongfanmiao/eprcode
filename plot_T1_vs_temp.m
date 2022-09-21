% plot_T1_vs_temp   plot T1 value against temperature
%
%   plot_T1_vs_temp(FitResult)
%   plot_T1_vs_temp(FitResult,'inverseT1')
%   plot_T1_vs_temp(FitResult,'Color','k','Marker','*','MarkerSize',50)
%
%   Input:
%       FitResult       table of fitting result
%       inverseT1       plot the inverse of T1 values
%       Color           color, 'blue' by default
%       Marker          marker, 'o' by default
%       MarkerSize      marker size, 100 by default
%       filled          filled circle or not
%       
function plot_T1_vs_temp(FitResult,varargin)

par = inputParser;
par.addParameter('Inverse', false, @islogical);
par.addParameter('LogScale', false, @islogical);
par.addParameter('Color', 'b', @(x)isstring(x)|ischar(x));
par.addParameter('Marker', 'o', @(x)any(validatestring(x, ...
    {'o','+','*','.','x','-','|','s','d','^','v','>','<','p','h'})));
par.addParameter('MarkerSize', 100, @isnumeric);
par.addParameter('filled', true, @islogical);
par.KeepUnmatched = true;

parse(par, varargin{:});
args = par.Results;

temp = FitResult.Temperature;

try
    T1 = FitResult.T1; % mono or stretch exponential
catch
    try
        T1 = FitResult.T_long; % bi-exponential
    catch
        error('There is no T1 or T_long in the FitResult table');
    end
end

if args.Inverse
    y = 1./T1;
    ylab = sprintf('{T_1}^{-1} (%s^{-1})', ...
        FitResult.Properties.VariableUnits{2});
    ttl = '{T_1}^{-1} versus Temperature';
else
    y = T1;
    ylab = sprintf('T_1 (%s)', ...
        FitResult.Properties.VariableUnits{2});
    ttl = 'T_1 versus Temperature';
end

additionalArgs = {};
if args.filled
    additionalArgs = {'filled'};
end
scatter(temp, y, args.MarkerSize, args.Marker, args.Color, ...
    additionalArgs{:});
xlabel('Temperature (K)');
ylabel(ylab);
title(ttl)
set(gcf,'color','w');
box on

if args.LogScale
    set(gca, 'XScale', 'log', 'YScale', 'log');

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end