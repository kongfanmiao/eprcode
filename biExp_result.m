function biExp_result(FitResult, plot_percentage)

% Plot the T_long, T_short and their percentages for varying temperatures
% FitResult: output of fit_T1(fit_Tm) with the bi-exponential fitting model

arguments
    FitResult   table
    plot_percentage logical = true
end


if plot_percentage
    yyaxis left
end

scatter(FitResult.Temperature, FitResult.T_long, 100, 'blue', 'filled','o');
hold on
scatter(FitResult.Temperature, FitResult.T_short, 100,'blue', 'o');
ylabel(sprintf("T_{long} & T_{short} (%s)", ...
    FitResult.Properties.VariableUnits{2}));

legend('T_{long}', 'T_{short}', "Location","best");

if plot_percentage
    yyaxis right
    scatter(FitResult.Temperature, ...
        FitResult.('T_long Percentage'), 100, 'filled', '^');
    scatter(FitResult.Temperature, ...
        FitResult.('T_short Percentage'), 100, 'v');
    ylabel('Percentage')
    
    legend('T_{long}', 'T_{short}', 'T_{long} percentage', ...
        'T_{short} percentage', "Location","best");
end

xlabel('Temperature (K)');
hold off
set(gcf,'color','w');
box on

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end
