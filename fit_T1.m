% fit_T1    Fit the T1 Picket Fence measurement data
%
%   FitResult = fit_T1(path, keywords)
%   FitResult = fit_T1(path, keywords, FitFunc)
%   FitResult = fit_T1(path, keywords, 'MarkerSize', 20)
%   FitResult = fit_T1(path, keywords, 'TimeUnit', 'ns')
%
%   Input:
%       path        directory to find files
%       keywords    keywords based on which to find files
%                   - string: "80K" (if there's only one keyword)
%                   - character vector: '80K' (if there's only one keyword)
%                   - string array: ["80K", "3480G"]
%                   - cell array of character vectors: {'80K', '3480G'}
%       FitFunc     fit function, usually use bi-exponential  
%                   - {'Exp1', 'Exp2', 'StrExp'}
%       MarkerSize  marker size in scatter plot, 20 by default
%                   - Name-value pair arguments
%       TimeUnit    X axis unit, use microsecond (us) by default 
%                   - {'ns', 'us', 'ms', 's'}
%       LegendColumns
%                   number of columns of legends
%       LegendPosition
%                   position of legend, eg. best, bestoutside
%       Color       use color map (default) or set color
%       Exclude     exclude certain temperatures
%       SeperatePlot
%                   plot each temperature in seperate panels
%       PlotTwoExp  plot the two processes in bi-exponential model
%       rcRange     range of rc in bi-exp fitting
%       rkRange     range of rk in bi-exp fitting
%
%   Output:
%       FitResult   table containing fitting results.
%
%   Example:
%       Res = fit_T1('./Data',["3480G","5K"],'Exp2','TimeUnit','us')
%


function FitResult = fit_T1(path, keywords, varargin)

% parse the input arguments
par = inputParser;

checkFitFunc = @(x) any(validatestring(x, {'Exp1', 'Exp2', 'StrExp'}));
par.addOptional('FitFunc', 'Exp2', checkFitFunc);

par.addParameter("MarkerSize", 10, @isnumeric);

checkTimeUnit = @(x) any(validatestring(x, {'ns', 'us', 'ms', 's'}));
par.addParameter("TimeUnit", 'ns', checkTimeUnit);
par.addParameter('LegendColumns', nan, @isnumeric);
par.addParameter('LegendPosition', 'southeast', @ischar);
par.addParameter('Color', nan, @ischar);
par.addParameter('Exclude', nan, @isnumeric);
par.addParameter('SeperatePlot', false, @islogical);
par.addParameter('PlotTwoExp', false, @islogical);
par.addParameter('rcRange', [0.01 100]);
% par.addParameter('rkRange', [0.01 100]);

par.KeepUnmatched = true;

parse(par, varargin{:});
args = par.Results;
fitfunc = args.FitFunc;

rc0 = args.rcRange(1);
rc1 = args.rcRange(2);
% rk0 = args.rkRange(1);
% rk1 = args.rkRange(2);

if args.PlotTwoExp
    % verify 1) SeperatePlot is true, 2) FitFunc is Exp2
    assert(args.SeperatePlot, ['SepeartePlot must be true in order to ' ...
        'plot two processes seperately in bi-exponential model']);
    assert(strcmp(fitfunc,'Exp2'), ['You must use bi-exponential model if ' ...
        'PlotTwoExp is true']);
end

fitfuncFullName = containers.Map( ...
    ["Exp1", "Exp2", "StrExp"], ...
    ["Mono-exponential", "Bi-exponential", "Stretched Exponential"]);

% Get the list of file names
% Convert to string if it's cell array of character vectors
keywords = string(keywords); 
keywords(end+1) = "T1PF";
keywords = unique(keywords);
filesTable = sort_temperature(path, keywords);
if ~isnan(args.Exclude)
    for i = 1:numel(args.Exclude)
        t = args.Exclude(i);
        filesTable(filesTable.Temperature==t,:) = [];
    end
end
Temperature = filesTable.Temperature;
Files = filesTable.Files;
numFiles = length(Files);

% Create the figure title based on keywords
titleStr = join(keywords, " ");
% Create figure legends
labels = cell(2*numFiles,1);

timeUnit = args.TimeUnit;
switch timeUnit
    case 'ns'
        timeScaleFactor = 1;
    case 'us'
        timeScaleFactor = 1e3;
        timeUnit = '\mus';
    case 'ms'
        timeScaleFactor = 1e6;
    case 's'
        timeScaleFactor = 1e9;
end

% Set the x and y axis label
xlabelStr = sprintf("Recovery time (%s)", timeUnit);
ylabelStr = "Signal (arb. u.)";

xMax = 0;
xMin = inf;
yMax = 0;
yMin = inf;
colorList = cool(numFiles);

% Initialize FitResult table
switch fitfunc
    case "Exp1"
        FitResult = table('VariableNames',{'Temperature', 'T1', 'T1ConfInt_lower', 'T1ConfInt_upper'}, ...
            'Size', [numFiles, 4], 'VariableTypes',repelem({'double'},4));
        FitResult.Properties.VariableUnits = {'K', timeUnit,'',''};
    case "Exp2"
        FitResult = table('VariableNames',{'Temperature', ...
            'T_long', 'T_long Weight', 'TlongConfInt_lower', 'TlongConfInt_upper', 'T_short', 'T_short Weight', 'TshortConfInt_lower', 'TshortConfInt_upper'}, ...
            'Size', [numFiles, 9], 'VariableTypes',repelem({'double'},9));
        FitResult.Properties.VariableUnits = {'K', timeUnit, '','','', ...
            timeUnit, '','',''};
    case "StrExp"
        FitResult = table('VariableNames',{'Temperature', ...
            'T1', 'Stretch Factor', 'T1ConfInt_lower', 'T1ConfInt_upper'}, ...
            'Size', [numFiles, 5], 'VariableTypes',repelem({'double'}, 5));
        FitResult.Properties.VariableUnits = {'K', timeUnit, '','',''};
end

FitResult.Temperature = Temperature;
FitResult.("R-Squared") = zeros(numFiles,1);

% Do the fitting
for i = 1:numFiles
    f = Files{i};
    [x, y] = eprload(fullfile(path, f));
    y = real(y);
    % change x axis for T1PF measurement
    x = change_x_axis(x, fullfile(path,f));
    % Change x axis unit
    x = x/timeScaleFactor;
    % Rescale y to be in the same scale as x
    y = y*max(x)/max(y);
    % first use exponfit to get an estimate
    % "c1 + c2*exp(-k*x)"
    [k,c,yfit] = exponfit(x,y,1);
    switch fitfunc
        case "Exp1"
            func = @(c1, c2, k, x) c1 + c2*exp(-x/k);
            fo = fitoptions(func,"StartPoint",[c 1/k]);
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            FitResult.T1(i) = curve.k;
            ci = confint(curve);
            FitResult.T1ConfInt_lower(i) = ci(1,3);
            FitResult.T1ConfInt_upper(i) = ci(2,3);
            FitResult.("R-Squared")(i) = gof.rsquare;
            yfit = curve(x);
        case "Exp2"
            func = @(c1, c2, rc, k1, k2, x) ...
                c1 + c2*exp(-x/k1) + c2*rc*exp(-x/(k2));
            % limiting slow process having larger weight
            bounds = [c(1)/2 c(2)/100 rc0 0.001/k 0.001/k;...
                      c(1)*2 c(2)*100 rc1 1000/k 1000/k];
            fo = fitoptions(func, ...
                "StartPoint",[c(1) c(2) 1 1/k 1.1/k], ...
                "Lower",min(bounds,[],1), ...
                "Upper",max(bounds,[],1));
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            if curve.k1 > curve.k2 % k1 is T long (Tl)
                TlPerc = abs(curve.c2)./(abs(curve.c2)+abs(curve.c2*curve.rc));
                TsPerc = abs(curve.c2*curve.rc)./(abs(curve.c2)+abs(curve.c2*curve.rc));
                ci = confint(curve);
                TlConfInt_lower = ci(1,4);
                TlConfInt_upper = ci(2,4);
                TsConfInt_lower = ci(1,5);
                TsConfInt_upper = ci(2,5);
                FitResult{i,2:9} = [curve.k1, TlPerc, TlConfInt_lower, TlConfInt_upper, ...
                    curve.k2, TsPerc, TsConfInt_lower, TsConfInt_upper];
            else % k2 is T long (Tl)
                TsPerc = abs(curve.c2)./(abs(curve.c2)+abs(curve.c2*curve.rc));
                TlPerc = abs(curve.c2*curve.rc)./(abs(curve.c2)+abs(curve.c2*curve.rc));
                ci = confint(curve);
                TsConfInt_lower = ci(1,4);
                TsConfInt_upper = ci(2,4);
                TlConfInt_lower = ci(1,5);
                TlConfInt_upper = ci(2,5);
                FitResult{i,2:9} = [curve.k2, TlPerc, TlConfInt_lower, TlConfInt_upper, ...
                    curve.k1, TsPerc, TsConfInt_lower, TsConfInt_upper];
            end
            FitResult.("R-Squared")(i) = gof.rsquare;
            yfit = curve(x);
        case "StrExp"
            func = "c1 + c2*exp(-(x/k)^c3)";
            bounds = [c/100 0.2 k/100;...
                      c*100 2.0 k*100]
            fo = fitoptions(func, ...
                "StartPoint",[c 1 1/k], ...
                "Lower",min(bounds,[],1), ...
                "Upper",max(bounds,[],1));
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            FitResult{i,2:3} = [curve.k, curve.c3];
            ci = confint(curve);
            FitResult.T1ConfInt_lower(i) = ci(1,4);
            FitResult.T1ConfInt_upper(i) = ci(2,4);
            FitResult.("R-Squared")(i) = gof.rsquare; 
            yfit = curve(x);
    end
    % Normalize y axis
    yM = max(y);
    y = y/yM;
    yfit = yfit/yM;
    if isnan(args.Color)
        c = colorList(i,:);
    else
        c = args.Color;
    end
    if ~args.SeperatePlot
        scatter(x,y,args.MarkerSize, "MarkerEdgeColor",c);
        labels{2*i-1} = strcat(num2str(Temperature(i))," K data");
        hold on
        plot(x,yfit, "Color", c);
        labels{2*i} = strcat(num2str(Temperature(i))," K Fit");
        hold on
    else
        figure();
        scatter(x,y, 30, "MarkerEdgeColor",'b');
        hold on
        plot(x,yfit, "Color", 'k', 'LineWidth',1);
        hold off
        box on
%         set(gca,'xscale','log');
        ylabel(ylabelStr);
        xlabel(xlabelStr);
        legend(strcat(num2str(Temperature(i))," K data"), ...
            strcat(num2str(Temperature(i))," K Fit"), ...
            'Location', 'best');
        title(sprintf('%s K', num2str(Temperature(i))));
        xlim([-0.01*max(x),1.01*max(x)]);
        yticks([]);
        if args.PlotTwoExp
            if T_long == curve.k1
                coef_long = curve.c2;
                coef_short = curve.c2*curve.rc;
            else
                coef_long = curve.c2*curve.rc;
                coef_short = curve.c2;
            end
            onlyslow = curve.c1 + coef_long*exp(-x/T_long);
            onlyfast = curve.c1 + coef_short*exp(-x/T_short);
            hold on
            plot(x, onlyslow/yM, 'Color', 'red', 'LineWidth',1);
            plot(x, onlyfast/yM, 'Color', 'green', 'LineWidth',1);
            legend(strcat(num2str(Temperature(i))," K data"), ...
                strcat(num2str(Temperature(i))," K Fit"), ... 
                sprintf(['only slow process (coef: %0.2e, T_{slow}: ' ...
                '%.2f %s)'], coef_long, T_long, timeUnit), ...
                sprintf(['only fast process (coef: %0.2e, T_{fast}: ' ...
                '%.2f %s)'], coef_short, T_short, timeUnit), ...
                'location','best');
        end
    end
    xMin = min(xMin, min(x));
    xMax = max(xMax, max(x));
    yMin = min(yMin, min(real(y)));
    yMax = max(yMax, max(real(y)));
end

if isnan(args.LegendColumns)
    numCol = ceil(numel(labels)/16);
else
    numCol = args.LegendColumns;
end

if ~args.SeperatePlot
    set(gca,'xscale','log');
    title(strcat(titleStr, sprintf("\n(%s)",fitfuncFullName(fitfunc))), ...
        "Interpreter","none");
    legend(labels, "Location", args.LegendPosition, "Interpreter","none", ...
        'NumColumns', numCol, 'FontSize',6);
    yticks([]);
    xlim([xMin xMax]);
    ylim([yMin yMax]);
    ylabel(ylabelStr);
    xlabel(xlabelStr);
    hold off
    box on
end
set(gcf,'color','w');
% axP = get(gca,'Position');  
% set(gca, 'Position', axP);

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end
