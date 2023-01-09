function biExp_result(FitResult, plot_percentage, xscale)

% Plot the T_long, T_short and their percentages for varying temperatures
% FitResult: output of fit_T1(fit_Tm) with the bi-exponential fitting model

arguments
    FitResult   table
    plot_percentage logical = true
    xscale char = 'linear'
end


if plot_percentage
    yyaxis left
end

scatter(FitResult.Temperature, FitResult.T_long, 100, 'blue', 'filled','o');
hold on
scatter(FitResult.Temperature, FitResult.T_short, 100, 'blue', 'o');
ylabel(sprintf("T_{long} & T_{short} (%s)", ...
    FitResult.Properties.VariableUnits{2}));
set(gca, 'ycolor', 'blue')
legend('T_{long}', 'T_{short}', "Location","best");

if plot_percentage
    yyaxis right
    scatter(FitResult.Temperature, ...
        FitResult.('T_long Weight'), 100, 'filled', '^');
    scatter(FitResult.Temperature, ...
        FitResult.('T_short Weight'), 100, 'v');
    ylabel('Weight')
    
    legend('T_{long}', 'T_{short}', 'T_{long} weight', ...
        'T_{short} weight', "Location","best");
end

xlabel('Temperature (K)');
hold off
set(gcf,'color','w');
set(gca, 'xscale', xscale);
box on

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end
