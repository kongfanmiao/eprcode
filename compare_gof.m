% compare_gof       compare goodness of fitting
%
%   compare_gof(FitResult1, label1, FitResult2, label2, ...)
%

function compare_gof(varargin)

if mod(nargin,2) == 1
    error("The arguments should be in pairs. Each FitResult table is " + ...
        "accompanied by a label. So you must provide even number of arguments");
end
FitResults = {varargin{1:2:end}};
labels = {varargin{2:2:end}};
temp = cellfun(@(t)t.Temperature, FitResults, 'UniformOutput',false);
rsquare = cellfun(@(t)t.("R-Squared"), FitResults, 'UniformOutput',false);

clf;
set(gcf,'color','w');
for i = 1:numel(temp)
    scatter(temp{i}, rsquare{i}, 'filled')
    hold on
end
legend(labels);
xlabel('Temperature (K)');
ylabel('R-Squared')
hold off
box on


% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end