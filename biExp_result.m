function biExp_result(FitResult, varargin)

% Plot the T_long, T_short and their percentages for varying temperatures
% FitResult: output of fit_T1(fit_Tm) with the bi-exponential fitting model

par = inputParser;
par.addParameter('PlotWeight', true, @islogical);
par.addParameter('XScale', 'linear', @ischar);
par.addParameter('YScale', 'linear', @ischar);
par.addParameter('SeperatePlot', true, @islogical);

par.KeepUnmatched = true;
parse(par, varargin{:});
args = par.Results;


if ~args.SeperatePlot % plot T1/Tm and weight in one graph

if args.PlotWeight % plot weight for fast and slow process
    yyaxis left
end

scatter(FitResult.Temperature, FitResult.T_long, 100, 'blue', 'filled','o');
hold on
scatter(FitResult.Temperature, FitResult.T_short, 100, 'blue', 'o');
ylabel(sprintf("T_{long} & T_{short} (%s)", ...
    FitResult.Properties.VariableUnits{2}));
set(gca, 'ycolor', 'blue', 'Yscale', args.YScale)
legend('T_{long}', 'T_{short}', "Location","best");

if args.PlotWeight
    yyaxis right
    scatter(FitResult.Temperature, ...
        FitResult.('T_long Weight'), 100, 'filled', '^');
    scatter(FitResult.Temperature, ...
        FitResult.('T_short Weight'), 100, 'v');
    ylabel('Weight')
    
    legend('T_{long}', 'T_{short}', 'T_{long} weight', ...
        'T_{short} weight', "Location","best");
end

else % plot T1/Tm and weight in two seperate panel

t = tiledlayout(2,1,'TileSpacing','Compact');

nexttile
scatter(FitResult.Temperature, FitResult.T_long, 40, ...
    'markerfacecolor','b', 'markeredgecolor','none');
hold on
scatter(FitResult.Temperature, FitResult.T_short, 40, ...
    'markerfacecolor','none', 'markeredgecolor',"#D95319");
ylabel(sprintf("T_{long} & T_{short} (%s)", ...
    FitResult.Properties.VariableUnits{2}));
set(gca,'XTick',[], 'Yscale', args.YScale);
legend('T_{long}', 'T_{short}', "Location","best");
box on

nexttile
scatter(FitResult.Temperature, ...
    FitResult.('T_long Weight'), 40, '^', ...
    'markerfacecolor','b', 'markeredgecolor','none');
hold on
scatter(FitResult.Temperature, ...
    FitResult.('T_short Weight'), 40, '^', ...
    'markerfacecolor','none', 'markeredgecolor',"#D95319");
ylabel('Weight');
legend('T_{long} weight', 'T_{short} weight', 'location', 'best');
ylim([0,1]);

end

xlabel('Temperature (K)');
hold off
set(gcf,'color','w');
set(gca, 'xscale', args.XScale);
box on

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end
