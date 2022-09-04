function plot_TmStretch(FitResult)

% Plot the Tm, and stretch factor for varying temperatures
% FitResult: output of fit_Tm with the stretched exponential fitting model

xlabel('Temperature (K)');

yyaxis left
scatter(FitResult.Temperature, FitResult.Tm, 100, 'filled','o','b');
ylabel(sprintf("T_m (%s)", FitResult.Properties.VariableUnits{2}));

yyaxis right
sf = FitResult.('Stretch Factor');
scatter(FitResult.Temperature, sf, 100, '*');
ylabel('Stretch Factor')
ylim([2*min(sf)-max(sf), 2*max(sf)-min(sf)]);

legend('T_m', 'Stretch Factor', "Location","best");
hold off
set(gcf,'color','w');
box on

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end
