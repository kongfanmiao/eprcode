function strExp_result(FitResult, xscale)

% Plot the Tm or T1, and stretch factor for varying temperatures
% FitResult: output of fit_Tm o rfit_T1 with the stretched exponential fitting model

arguments
    FitResult table
    xscale char = 'linear'
end

xlabel('Temperature (K)');

getVarName = @(x) inputname(1);
try
    T = FitResult.Tm;
    label = 'T_m';
catch
    try
        T = FitResult.T1;
        label = 'T_1';
    catch
        error(sprintf('%s must contains either Tm or T1', ...
            getVarName(FitResult)));
    end
end
yyaxis left
scatter(FitResult.Temperature, T, 100, 'filled','o','b');
ylabel(sprintf("%s (%s)", label, FitResult.Properties.VariableUnits{2}));
ylim([min(T)*0.9 max(T)*1.1])
xlim([min(FitResult.Temperature)*0.9 max(FitResult.Temperature)*1.1])
set(gca, 'ycolor', 'blue')
yyaxis right
sf = FitResult.('Stretch Factor');
scatter(FitResult.Temperature, sf, 100, '*');
ylabel('Stretch Factor');
ylim([2*min(sf)-max(sf), 2*max(sf)-min(sf)]);

legend(label, 'Stretch Factor', "Location","best");
hold off
set(gcf,'color','w');
set(gca, 'xscale', xscale);
box on
grid on
% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end
