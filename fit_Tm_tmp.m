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
par.addParameter('FitModulation', false, @islogical);
par.addParameter('TwoOmega', false, @islogical);

par.addParameter('MarkerSize', 20, @isnumeric);
par.addParameter('DataColor', nan, @(x) ischar(x)|isstring(x));
par.addParameter('FitColor', nan, @(x) ischar(x)|isstring(x));

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
        if ~args.FitModulation
            FitResult = table( ...
                'VariableNames',{'Temperature', 'Tm', 'R-Squared'}, ...
                'Size', [numFiles, 3], ...
                'VariableTypes',repelem({'double'},3));
            FitResult.Properties.VariableUnits = {'K', timeUnit, ''};
        else
            if ~args.TwoOmega
                FitResult = table( ...
                    'VariableNames',{'Temperature', 'Tm', 'C_mod', ...
                        'Modulation Frequency', 'R-Squared'}, ...
                    'Size', [numFiles, 5], ...
                    'VariableTypes',repelem({'double'},5));
                FitResult.Properties.VariableUnits = {
                    'K', timeUnit, '', 'MHz', ''};
            else
                FitResult = table( ...
                    'VariableNames', ...
                        {'Temperature', 'Tm', 'C_mod1', 'C_mod2', ...
                        'Modulation Frequency', 'R-Squared'}, ...
                    'Size', [numFiles, 6], ...
                    'VariableTypes',repelem({'double'},6));
                FitResult.Properties.VariableUnits = {
                    'K', timeUnit, '', '', 'MHz', ''};
            end
        end
    case "Exp2"
        if ~args.FitModulation
            FitResult = table('VariableNames', ...
                {'Temperature', 'T_long', 'T_long Percentage', ...
                'T_short', 'T_short Percentage', 'R-squared'}, ...
                'Size', [numFiles, 6], ...
                'VariableTypes',repelem({'double'},6));
            FitResult.Properties.VariableUnits = {
                'K', timeUnit, '', timeUnit, '', ''};
        else
            if ~args.TwoOmega
                FitResult = table('VariableNames', ...
                    {'Temperature', 'T_long', 'T_long Percentage', ...
                        'T_short', 'T_short Percentage', 'C_mod', ...
                        'Modulation Frequency', 'R-squared'}, ...
                    'Size', [numFiles, 8], ...
                    'VariableTypes',repelem({'double'},8));
                FitResult.Properties.VariableUnits = {
                    'K', timeUnit, '', timeUnit, '', '', 'MHz', ''};
            else
                FitResult = table('VariableNames', ...
                    {'Temperature', 'T_long', 'T_long Percentage', ...
                        'T_short', 'T_short Percentage', 'C_mod1', ...
                        'C_mod2','Modulation Frequency', 'R-squared'}, ...
                    'Size', [numFiles, 9], ...
                    'VariableTypes',repelem({'double'},9));
                FitResult.Properties.VariableUnits = {
                    'K', timeUnit, '', timeUnit, '', '', '', 'MHz', ''};
            end
        end
    case "StrExp"
        if ~args.FitModulation
            FitResult = table('VariableNames', ...
                {'Temperature', 'Tm', 'Stretch Factor', 'R-squared'}, ...
            'Size', [numFiles, 4], ...
            'VariableTypes',repelem({'double'},4));
            FitResult.Properties.VariableUnits = {'K', timeUnit, '', ''};
        else
            if ~args.TwoOmega
                FitResult = table('VariableNames', ...
                    {'Temperature', 'Tm', 'Stretch Factor', ...
                    'C_mod', 'Modulation Frequency', 'R-Squared'}, ...
                    'Size', [numFiles, 6], ...
                    'VariableTypes',repelem({'double'},6));
                FitResult.Properties.VariableUnits = ...
                    {'K', timeUnit, '', '', 'MHz', ''};
            else
                FitResult = table( ...
                    'VariableNames', ...
                        {'Temperature', 'Tm', 'Stretch Factor', 'C_mod1', ...
                        'C_mod2', 'Modulation Frequency', 'R-Squared'}, ...
                    'Size', [numFiles, 7], ...
                    'VariableTypes',repelem({'double'},7));
                FitResult.Properties.VariableUnits = ...
                    {'K', timeUnit, '', '', '', 'MHz', ''};
            end
        end
end

    switch fitfunc
        case "Exp1"
            if ~args.FitModulation
                func = @(c1, c2, k, x) c1 + c2*exp(-x/k);
                bounds = [c(1)/2,c(2)/100, 0.5/k;
                          c(1)*2,c(2)*100, 2/k];
                startPoint = [c(1), c(2), 1/k];
            else
                if ~args.TwoOmega
                    func = @(c1, c2, k, c3, omega, phi_1, x) ...
                        (c1 + c2*exp(-x/k)).*(1 + c3*sin(omega*x+phi_1));
                    bounds = [c(1)/2,c(2)/100, 0.5/k,1e-3,omega/10,-inf;
                              c(1)*2,c(2)*100, 2/k,100,omega*10,inf];
                    startPoint = [c(1), c(2), 1/k, 1, omega, 0];
                else
                    func = @(c1, c2, k, c3, c4, omega, phi_1, phi_2, x) ...
                        (c1 + c2*exp(-x/k)).*(1 + c3*sin(omega*x+phi_1) + ...
                        c4*sin(2*omega*x+phi_2));
                    bounds = [c(1)/2,c(2)/100, 0.5/k,1e-3,1e-3, ...
                                                    omega/10,-inf,-inf;
                              c(1)*2,c(2)*100, 2/k,100,100, ...
                                                    omega*10,inf,inf];
                    startPoint = [c(1), c(2), 1/k, 1, 1, omega, 0, 0];
                end
            end


            fo = fitoptions(func, ...
                "StartPoint",startPoint, ...
                "Lower",min(bounds,[],1), ...
                "Upper",max(bounds,[],1));
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            FitResult.Tm(i) = curve.k;
            FitResult.("R-Squared")(i) = gof.rsquare;
            yfit = curve(x);
        case "Exp2"
            func = @(c1, c2, c3, k1, k2, x) ...
                c1 + c2*exp(-x/k1) + c3*exp(-x/k2);
            bounds = [c(1)/2 c(2)/100 c(2)/100 0.1/k 0.0001/k;...
                      c(1)*2 c(2)*100 c(2)*100 1000/k 10/k];
            fo = fitoptions(func, ...
                "StartPoint",[c(1) c(2) c(2) 1/k 1/k], ...
                "Lower",min(bounds,[],1), ...
                "Upper",max(bounds,[],1));
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            [T_long, T_longPerc, T_short, T_shortPerc] = ...
                get_Tlong_and_Tshort([curve.k1, curve.k2], ...
                [curve.c1, curve.c2, curve.c3]);
            FitResult{i,2:end-1} = [T_long, T_longPerc, T_short, T_shortPerc];
            FitResult.("R-Squared")(i) = gof.rsquare; 
            yfit = curve(x);
        case "StrExp"
            func = @(c1, c2, c3, k, x) ...
                c1 + c2*exp(-(x/k).^c3);
            bounds = [c/100 0.7 0.01/k;...
                      c*100 1.3 100/k];
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

FitResult.Temperature = Temperature;

% Do the fitting
    function [curve, gof] = fit_single_trace(f,fitfunc)
        [x, y] = eprload(fullfile(path, f));
        % Change x axis unit
        x = x/timeScaleFactor;
        y = real(y);          
        % first use exponfit to get an estimate
        % "c1 + c2*exp(-k*x)"
        [k,c,yfit] = exponfit(x,y,1);
        if args.FitModulation
            % use fft to get an estimate of the modulation frequency
            N = numel(x);
            yFFT = abs(fft(y-yfit))/N;
            freqAxis = (0:N/2-1)/(N*(x(2)-x(1)));
            yFFT = yFFT(1:N/2);
            [~, loc] = max(yFFT);
            modFreq = freqAxis(loc); % Depends on time unit
            omega = modFreq*2*pi;
        end
        fo = fitoptions(fitfunc, ...
            "StartPoint",[c 1 1/k], ...
            "Lower",min(bounds,[],1), ...
            "Upper",max(bounds,[],1));
        ft = fittype(fitfunc,"options",fo);
        [curve, gof] = fit(x,y,ft); 
        yfit = curve(x);

        % Normalize y axis
        yM = max(y);
        y = y/yM;
        yfit = yfit/yM;
    end

    if ~isnan(args.DataColor)
        dataC = args.DataColor;
    else
        dataC = colorList(i,:);
    end
    scatter(x,y,args.MarkerSize, "MarkerEdgeColor", dataC);
    labels{2*i-1} = strcat(num2str(Temperature(i))," K data");
    hold on
    if ~isnan(args.FitColor)
        fitC = args.FitColor;
    else
        fitC = colorList(i,:);
    end
    plot(x,yfit, "Color", fitC);
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

