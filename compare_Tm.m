% compare_Tm       compare goodness of fitting
%
%   compare_Tm(FitResult1, label1, FitResult2, label2, ...)
%

function compare_Tm(varargin)

if mod(nargin,2) == 1
    error("The arguments should be in pairs. Each FitResult table is " + ...
        "accompanied by a label. So you must provide even number of arguments");
end
FitResults = {varargin{1:2:end}};
labels = {varargin{2:2:end}};
temp = cellfun(@(t)t.Temperature, FitResults, 'UniformOutput',false);
Tmall = {};
for i = 1:numel(FitResults)
    FitResult = FitResults{i};
    try
        Tm = FitResult.Tm; % mono or stretch exponential
    catch
        try
            Tm = FitResult.T_long; % bi-exponential
        catch
            error('There is no Tm or T_long in the FitResult table');
        end
    end
    Tmall{end+1} = Tm;
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
    scatter(temp{i}, Tmall{i}, 50, marker, color)
    hold on
end
legend(labels);
xlabel('Temperature (K)');
ylabel(sprintf('Tm (%s)', FitResult.Properties.VariableUnits{2}));
title('Compare Tm from different fitting model')
hold off
box on


% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end