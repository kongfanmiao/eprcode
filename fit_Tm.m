% fit_Tm    Fit the Tm measurement data
%
%   FitResult = fit_Tm(path, keywords)
%   FitResult = fit_Tm(path, keywords, FitFunc)
%   FitResult = fit_Tm(path, keywords, 'MarkerSize', 20)
%   FitResult = fit_Tm(path, keywords, 'TimeUnit', 'ns')
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
%
%   Output:
%       FitResult   table containing fitting results.
%
%   Example:
%       Res = fit_T1('./Data',["3480G","5K"],'StrExp','TimeUnit','ns')
%


function FitResult = fit_Tm(path, keywords, varargin)

% parse the input arguments
par = inputParser;

checkFitFunc = @(x) any(validatestring(x, {'Exp1', 'Exp2', 'StrExp'}));
par.addOptional('FitFunc', 'Exp2', checkFitFunc);

par.addParameter('MarkerSize', 20, @isnumeric);

checkTimeUnit = @(x) any(validatestring(x, {'ns', 'us', 'ms', 's'}));
par.addParameter('TimeUnit', 'ns', checkTimeUnit);

par.KeepUnmatched = true;

parse(par, varargin{:});
args = par.Results;
fitfunc = args.FitFunc;

fitfuncFullName = containers.Map( ...
    ["Exp1", "Exp2", "StrExp"], ...
    ["Mono-exponential", "Bi-exponential", "Stretched Exponential"]);

% Get the list of file names
% Convert to string if it's cell array of character vectors
keywords = string(keywords); 
keywords(end+1) = "Tm";
keywords = unique(keywords);
filesTable = sort_temperature(path, keywords);
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
xlabelStr = sprintf(strcat("2",char([0xD835 0xDF0F]), " (%s)"), timeUnit);
ylabelStr = "Signal";

xMax = 0;
xMin = inf;
yMax = 0;
yMin = inf;
colorList = jet(numFiles);

% Initialize FitResult table
switch fitfunc
    case "Exp1"
        FitResult = table('VariableNames',{'Temperature', 'Tm'}, ...
            'Size', [numFiles, 2], 'VariableTypes',repelem({'double'},2));
        FitResult.Properties.VariableUnits = {'K', timeUnit};
    case "Exp2"
        FitResult = table('VariableNames',{'Temperature', ...
            'Tm', 'Tm Percentage', 'Tf', 'Tf Percentage'}, ...
            'Size', [numFiles, 5], 'VariableTypes',repelem({'double'},5));
        FitResult.Properties.VariableUnits = {'K', timeUnit, '', ...
            timeUnit, ''};
    case "StrExp"
        FitResult = table('VariableNames',{'Temperature', ...
            'Tm', 'Stretch Factor'}, ...
            'Size', [numFiles, 3], 'VariableTypes',repelem({'double'},3));
        FitResult.Properties.VariableUnits = {'K', timeUnit, ''};
end

FitResult.Temperature = Temperature;
FitResult.("R-Squared") = zeros(numFiles,1);

% Do the fitting
for i = 1:numFiles
    f = Files{i};
    [x, y] = eprload(fullfile(path, f));
    y = real(y);          
    % first use exponfit to get an estimate
    % "c1 + c2*exp(-k*x)"
    [k,c,yfit] = exponfit(x,y,1);
    switch fitfunc
        case "Exp1"
            func = "c1 + c2*exp(-x/k)";
            fo = fitoptions(func,"StartPoint",[c 1/k]);
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            FitResult.Tm(i) = curve.k;
            FitResult.("R-Square")(i) = gof.rsquare;
            yfit = curve(x);
        case "Exp2"
            func = "c1 + c2*exp(-x/k1) + c3*exp(-x/k2)";
            bounds = [c(1)/2 c(2)/100 c(2)/100 0.001/k 0.001/k;...
                      c(1)*2 c(2)*100 c(2)*100 1000/k 1000/k];
            fo = fitoptions(func, ...
                "StartPoint",[c(1) c(2) c(2) 1/k 1/k], ...
                "Lower",min(bounds,[],1), ...
                "Upper",max(bounds,[],1));
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            [Tm, TmPerc, Tf, TfPerc] = get_TandTf( ...
                [curve.k1, curve.k2], [curve.c1, curve.c2, curve.c3], ...
                timeScaleFactor);
            FitResult{i,2:end-1} = [Tm, TmPerc, Tf, TfPerc];
            FitResult.("R-Squared")(i) = gof.rsquare; 
            yfit = curve(x);
        case "StrExp"
            func = "c1 + c2*exp(-(x/k)^c3)";
            bounds = [c/100 0.2 0.01/k;...
                      c*100 2.0 100/k];
            fo = fitoptions(func, ...
                "StartPoint",[c 1 1/k], ...
                "Lower",min(bounds,[],1), ...
                "Upper",max(bounds,[],1));
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            FitResult{i,2:end-1} = [curve.k, curve.c3];
            FitResult.("R-Squared")(i) = gof.rsquare; 
            yfit = curve(x);
    end
    % Change x axis unit
    x = x/timeScaleFactor;
    % Normalize y axis
    yM = max(y);
    y = y/yM;
    yfit = yfit/yM;
    scatter(x,y,args.MarkerSize, "MarkerEdgeColor",colorList(i,:));
    labels{2*i-1} = strcat(num2str(Temperature(i))," K data");
    hold on
    plot(x,yfit, "Color", colorList(i,:));
    labels{2*i} = strcat(num2str(Temperature(i))," K Fit");
    hold on
    xMin = min(xMin, min(x));
    xMax = max(xMax, max(x));
    yMin = min(yMin, min(real(y)));
    yMax = max(yMax, max(real(y)));
end
title(strcat(titleStr, sprintf("\n(%s)",fitfuncFullName(fitfunc))), ...
    "Interpreter","none");
legend(labels, "Location","best", "Interpreter","none");
ylabel(ylabelStr);
xlabel(xlabelStr);
xlim([xMin xMax]);
ylim([yMin yMax]);
yticks([]);
hold off
set(gcf,'color','w');
box on

fig = gcf;
if fig.WindowStyle ~= "docked"
    set(fig,'position',[10,10,900,600]);
end

end

