function plot_simulation(path, fileName)

filePath = fullfile(path, fileName);

data = readtable(filePath, 'Range','A3');
names = readcell(filePath, 'Range', 'A1:C1');
units = readcell(filePath, 'Range', 'A2:C2');
data.Properties.VariableNames = names;
data.Properties.VariableUnits = units;

B = data{:,1};
spc = data{:,2};
fitSpc = data{:,3};

plot(B, spc, '-k',  'LineWidth', 1);
hold on
plot(B, fitSpc, '--b', 'LineWidth', 2);
xlim([min(B) max(B)]);
ylim([min(spc) max(max(spc), max(fitSpc))*1.1]);
xlabel(sprintf('%s (%s)', names{1}, units{1}));
ylabel(sprintf('%s (%s)', names{2}, units{2}));
titleStr = fileName(1:end-4);
title(titleStr, 'Interpreter','none');
legend(names{2:end});
yticks([]);
set(gca, 'LineWidth', 1.5);
set(gcf,'color','w');
box on

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end