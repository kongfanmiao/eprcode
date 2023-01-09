% compare_T1       compare goodness of fitting
%
%   compare_T1(FitResult1, label1, FitResult2, label2, ...)
%

function compare_T1(varargin)

if mod(nargin,2) == 1
    error("The arguments should be in pairs. Each FitResult table is " + ...
        "accompanied by a label. So you must provide even number of arguments");
end
FitResults = {varargin{1:2:end}};
labels = {varargin{2:2:end}};
temp = cellfun(@(t)t.Temperature, FitResults, 'UniformOutput',false);
T1all = {};
for i = 1:numel(FitResults)
    FitResult = FitResults{i};
    try
        T1 = FitResult.T1; % mono or stretch exponential
    catch
        try
            T1 = FitResult.T_long; % bi-exponential
        catch
            error('There is no T1 or T_long in the FitResult table');
        end
    end
    T1all{end+1} = T1;
end

clf;
set(gcf,'color','w');
for i = 1:numel(temp)
    label = labels{i};
    if contains(label, 'mono', 'IgnoreCase',true)
        color = 'r';
        marker = '*';
    elseif contains(label, 'bi', 'IgnoreCase',true)
        color = 'b';
        marker = '^';
    elseif contains(label, 'str', 'IgnoreCase',true)
        color = 'm';
        marker = 'v';
    else
        color = 'k';
        marker = 'o';
    end
    scatter(temp{i}, T1all{i}, 50, marker, color)
    hold on
end
legend(labels);
xlabel('Temperature (K)');
ylabel(sprintf('T1 (%s)', FitResult.Properties.VariableUnits{2}));
set(gca, 'yscale', 'log')
title('Compare T1 from different fitting model')
hold off
box on


% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end