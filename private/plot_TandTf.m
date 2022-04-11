function plot_TandTf(FitResult, label)

% Plot the T1(Tm), Tf and their percentages for varying temperatures
% FitResult: output of fit_T1(fit_Tm) with the bi-exponential fitting model

arguments
    FitResult   table
    label       char
end

switch label
    case 'T1'
        lb1 = 'T1';
        lb2 = 'T_1';
    case 'Tm'
        lb1 = 'Tm';
        lb2 = 'T_m';
end

xlabel('Temperature (K)');

yyaxis left
scatter(FitResult.Temperature, FitResult.(lb1), 100, 'filled','o');
hold on
scatter(FitResult.Temperature, FitResult.Tf, 100, 'o');
ylabel(sprintf("%s & T_f (%s)", ...
    lb2, FitResult.Properties.VariableUnits{2}));

yyaxis right
scatter(FitResult.Temperature, ...
    FitResult.(sprintf('%s Percentage', lb1)), 100, 'filled', '^');
scatter(FitResult.Temperature, FitResult.("Tf Percentage"), ...
    100, 'v');
ylabel('Percentage')

legend(lb2, 'T_f', sprintf('%s percentage', lb2), ...
    'T_f percentage', "Location","best");
hold off
set(gcf,'color','w');
box on

fig = gcf;
if fig.WindowStyle ~= "docked"
    set(fig,'position',[10,10,900,600]);
end

end
