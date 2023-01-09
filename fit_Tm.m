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
%       DataColor   data color
%       FitColor    fit color
%       TimeUnit    X axis unit, use microsecond (us) by default 
%                   - {'ns', 'us', 'ms', 's'}
%       StartingPoint
%       LowerBound
%       UpperBound
%       PrintFitCurve   true or flase
%       LegendColumns   number of columns
%       LegendPosition  'best'
%       yoffset     offset of each trace. 0 by default.
%       ModFreqRange  range of the modulation frequency
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

par.addParameter('MarkerSize', 10, @isnumeric);
par.addParameter('DataColor', nan, @(x) ischar(x)|isstring(x));
par.addParameter('FitColor', nan, @(x) ischar(x)|isstring(x));

checkTimeUnit = @(x) any(validatestring(x, {'ns', 'us', 'ms', 's'}));
par.addParameter('TimeUnit', 'ns', checkTimeUnit);

par.addParameter('StartingPoint', nan, @isnumeric);
par.addParameter('LowerBound', nan, @isnumeric);
par.addParameter('UpperBound', nan, @isnumeric);

par.addParameter('PrintFitCurve', false, @islogical);
par.addParameter('LegendColumns', nan, @isnumeric);
par.addParameter('LegendPosition', 'best', @ischar);
par.addParameter('yoffset', 0, @isnumeric);
par.addParameter('ModFreqRange',nan,@(x) isnumeric(x) && length(x)==2);

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
ylabelStr = "Signal (arb. u.)";

xMax = 0;
xMin = inf;
yMax = 0;
yMin = inf;
colorList = cool(numFiles);

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
                {'Temperature', 'T_long', 'T_long Weight', ...
                'T_short', 'T_short Weight', 'R-Squared'}, ...
                'Size', [numFiles, 6], ...
                'VariableTypes',repelem({'double'},6));
            FitResult.Properties.VariableUnits = {
                'K', timeUnit, '', timeUnit, '', ''};
        else
            if ~args.TwoOmega
                FitResult = table('VariableNames', ...
                    {'Temperature', 'T_long', 'T_long Weight', ...
                        'T_short', 'T_short Weight', 'C_mod', ...
                        'Modulation Frequency', 'R-Squared'}, ...
                    'Size', [numFiles, 8], ...
                    'VariableTypes',repelem({'double'},8));
                FitResult.Properties.VariableUnits = {
                    'K', timeUnit, '', timeUnit, '', '', 'MHz', ''};
            else
                FitResult = table('VariableNames', ...
                    {'Temperature', 'T_long', 'T_long Weight', ...
                        'T_short', 'T_short Weight', 'C_mod1', ...
                        'C_mod2','Modulation Frequency', 'R-Squared'}, ...
                    'Size', [numFiles, 9], ...
                    'VariableTypes',repelem({'double'},9));
                FitResult.Properties.VariableUnits = {
                    'K', timeUnit, '', timeUnit, '', '', '', 'MHz', ''};
            end
        end
    case "StrExp"
        if ~args.FitModulation
            FitResult = table('VariableNames', ...
                {'Temperature', 'Tm', 'Stretch Factor', 'R-Squared'}, ...
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

FitResult.Temperature = Temperature;

% Do the fitting
for i = 1:numFiles
    f = Files{i};
    [x, y] = eprload(fullfile(path, f));
    % Change x axis unit
    x = x/timeScaleFactor;
    y = real(y);      
    % first use exponfit to get an estimate
    % "c1 + c2*exp(-k*x)"
    [k,c,yfit] = exponfit(x,y,1);
    if args.FitModulation
        % use fft to get an estimate of the modulation frequency
        func = @(c1, c2, k, x) c1 + c2*exp(-x/k);
        bounds = [c(1)/2,c(2)/100, 0.5/k;
                  c(1)*2,c(2)*100, 2/k];
        fo = fitoptions(func, ...
            "StartPoint",[c(1), c(2), 1/k], ...
            "Lower",min(bounds,[],1), ...
            "Upper",max(bounds,[],2));
        ft = fittype(func,"options",fo);
        curve = fit(x,y,ft);
        yfit = curve(x);
        N = numel(x);
        yFFT = abs(fft(y-yfit))/N;
        freqAxis = (0:N/2-1)/(N*(x(2)-x(1)));
        yFFT = yFFT(1:N/2);
        [~, loc] = max(yFFT);
        modFreq = freqAxis(loc); % Depends on time unit
        omega = modFreq*2*pi;
        if ~isnan(args.ModFreqRange)
            omegaMin = args.ModFreqRange(1)*pi*timeScaleFactor/1e3;
            omegaMax = args.ModFreqRange(2)*pi*timeScaleFactor/1e3;
        else
            omegaMin = omega/10;
            omegaMax = omega*10;
        end
    end

    switch fitfunc
        case "Exp1"
            if ~args.FitModulation
                func = @(c1, c2, k, x) c1 + c2*exp(-x/k);
                bounds = [c(1)/2,c(2)/100, 0.5/k;
                          c(1)*2,c(2)*100, 2/k];
                startPoint = [c(1), c(2), 1/k];
                [curve, gof] = do_fit(x, y, func, startPoint, bounds);
            else
                if ~args.TwoOmega
                    func = @(c1, c2, k, c3, omega, phi, x) ...
                        c1 + c2*exp(-x/k).*(1 + c3*sin(omega*x+phi));
                    bounds = [c(1)/2,c(2)/100,0.5/k,1e-3,omegaMin,-inf;
                              c(1)*2,c(2)*100,2/k,100,omegaMax,inf];
                    startPoint = [c(1), c(2), 1/k, 1, omega, 0];
                    [curve, gof] = do_fit(x, y, func, startPoint, bounds);
                    FitResult.C_mod(i) = curve.c3;
                else
                    func = @(c1, c2, k, c3, c4, omega, phi_1, phi_2, x) ...
                        c1 + c2*exp(-x/k).*(1 + c3*sin(omega*x+phi_1) + ...
                        c4*sin(2*omega*x+phi_2));
                    bounds = [c(1)/2,c(2)/100, 0.1/k,1e-3,1e-3, ...
                                                    omegaMin,-inf,-inf;
                              c(1)*2,c(2)*100, 10/k,100,100, ...
                                                    omegaMax,inf,inf];
                    startPoint = [c(1), c(2), 1/k, 1, 1, omega, 0, 0];
                    [curve, gof] = do_fit(x, y, func, startPoint, bounds);
                    FitResult.C_mod1(i) = curve.c3;
                    FitResult.C_mod2(i) = curve.c4;
                end
                % the x here is 2tau. So omega*2 is the real frequency
                % so modFreq(real) = omega*2 /(2*pi)
                FitResult.('Modulation Frequency')(i) = ...
                    (curve.omega/pi)*1e3/timeScaleFactor; % MHz
            end
            FitResult.Tm(i) = curve.k;

        case "Exp2"
            if ~args.FitModulation
                func = @(c1, c2, c3, k1, r, x) ...
                        c1 + c2*exp(-x/k1) + c3*exp(-x/(k1*r));
                bounds = [c(1)/2,c(2)/100,c(2)/100,0.001/k,1;...
                          c(1)*2,c(2)*100,c(2)*100,1000/k,100];
                startPoint = [c(1), c(2), c(2), 1/k, 1];
                [curve, gof] = do_fit(x, y, func, startPoint, bounds);
            else
                if ~args.TwoOmega
                    func = @(c1, c2, c3, k1, r, c4, omega, phi, x) ...
                            c1 + (c2*exp(-x/k1) + c3*exp(-x/(k1*r))).*( ...
                            1 + c4*sin(omega*x+phi));
                    bounds = [c(1)/2,c(2)/100,c(2)/100,0.001/k,1, ...
                                                1e-3,omegaMin,-inf;
                              c(1)*2,c(2)*100,c(2)*100,1000/k,100, ...
                                                100,omegaMax,inf];
                    startPoint = [c(1), c(2), c(2), 1/k, 1, 1, omega, 0];
                    [curve, gof] = do_fit(x, y, func, startPoint, bounds);
                    FitResult.C_mod(i) = curve.c4;
                else
                    func = @(c1, c2, c3, k1, r, c4, c5, omega, phi_1, ...
                            phi_2, x) ...
                            c1 + (c2*exp(-x/k1) + c3*exp(-x/(k1*r))).*(1 + ...
                            c4*sin(omega*x+phi_1) + ...
                            c5*sin(2*omega*x+phi_2));
                    bounds = [c(1)/2,c(2)/100,c(2)/100,0.001/k,1, ...
                                            1e-3,1e-3,omegaMin,-inf,-inf;
                              c(1)*2,c(2)*100,c(2)*100,1000/k,100, ...
                                            100,100, omegaMax,inf,inf];
                    startPoint = [c(1), c(2), c(2), 1/k, 1, 1, 1, ...
                                            omega, 0, 0];
                    [curve, gof] = do_fit(x, y, func, startPoint, bounds);
                    FitResult.C_mod1(i) = curve.c4;
                    FitResult.C_mod2(i) = curve.c5;
                end
                FitResult.('Modulation Frequency')(i) = ...
                    (curve.omega/pi)*1e3/timeScaleFactor; % MHz
            end
            [T_long, T_longPerc, T_short, T_shortPerc] = ...
                get_Tlong_and_Tshort([curve.k1, curve.r*curve.k1], ...
                [curve.c1, curve.c2, curve.c3]);
            FitResult{i,2:5} = [T_long, T_longPerc, T_short, T_shortPerc];
        case "StrExp"
            if ~args.FitModulation
                func = @(c1, c2, beta, k, x) ...
                        c1 + c2*exp(-(x/k).^beta);
                bounds = [c(1)/2,c(2)/100,0.5,0.1/k;...
                          c(1)*2,c(2)*100,1.5,10/k];
                startPoint = [c(1),c(2),1,1/k];
                [curve, gof] = do_fit(x, y, func, startPoint, bounds);
            else
                if ~args.TwoOmega
                    func = @(c1, c2, beta, k, c3, omega, phi, x) ...
                            c1 + (c2*exp(-(x/k).^beta)).*(1 + ...
                            c3*sin(omega*x+phi));
                    bounds = [c(1)/2,c(2)/100,0.5,0.1/k,1e-3,omegaMin,-inf;
                              c(1)*2,c(2)*100,1.5,10/k,100,omegaMax,inf];
                    startPoint = [c(1),c(2),1,1/k,1,omega,0];
                    [curve, gof] = do_fit(x, y, func, startPoint, bounds);
                    FitResult.C_mod(i) = curve.c3;
                else
                    func = @(c1, c2, beta, k, c3, c4, omega, ...
                            phi_1, phi_2, x) ...
                            c1 + (c2*exp(-(x/k).^beta)).*(1 + ...
                            c3*sin(omega*x+phi_1) + c4*sin(2*omega*x+phi_2));
                    bounds = [c(1)/2,c(2)/100,0.5,0.1/k,1e-3,1e-3, ...
                                                    omegaMin,-inf,-inf;
                              c(1)*2,c(2)*100,1.5,10/k,100,100, ...
                                                    omegaMax,inf,inf];
                    startPoint = [c(1),c(2),1,1/k, 1, 1, omega, 0, 0];
                    [curve, gof] = do_fit(x, y, func, startPoint, bounds);
                    FitResult.C_mod1(i) = curve.c3;
                    FitResult.C_mod2(i) = curve.c4;
                end
                FitResult.('Modulation Frequency')(i) = ...
                    (curve.omega/pi)*1e3/timeScaleFactor; % MHz
            end
            FitResult.Tm(i) = curve.k;
            FitResult.('Stretch Factor')(i) = curve.beta;
    end
    if args.PrintFitCurve
        disp(curve);
    end
    yfit = curve(x);
    FitResult.("R-Squared")(i) = gof.rsquare;
    % Normalize label axis
    yM = max(y);
    y = y/yM;
    yfit = yfit/yM;
    if ~isnan(args.DataColor)
        dataC = args.DataColor;
    else
        dataC = colorList(i,:);
    end
    scatter(x,y+i*args.yoffset,args.MarkerSize, "MarkerEdgeColor", dataC);
    labels{2*i-1} = strcat(num2str(Temperature(i))," K data");
    hold on
    if ~isnan(args.FitColor)
        fitC = args.FitColor;
    else
        fitC = colorList(i,:);
    end
    plot(x,yfit+i*args.yoffset, "Color", fitC);
    labels{2*i} = strcat(num2str(Temperature(i))," K Fit");
    hold on
    xMin = min(xMin, min(x));
    xMax = max(xMax, max(x));
    yMin = min(yMin, min(real(y+i*args.yoffset)));
    yMax = max(yMax, max(real(y+i*args.yoffset)));

end
    
function [curve, gof] = do_fit(x,y,func,startPoint, bounds)
    if ~isnan(args.StartingPoint)
        sp = args.StartingPoint;
    else
        sp = startPoint;
    end
    if ~isnan(args.LowerBound)
        lb = args.LowerBound;
    else
        lb = min(bounds,[],1);
    end
    if ~isnan(args.UpperBound)
        ub = args.UpperBound;
    else
        ub = max(bounds,[],1);
    end
    fo = fitoptions(func, ...
        "StartPoint",sp, ...
        "Lower",lb, ...
        "Upper",ub);
    ft = fittype(func,"options",fo);
    [curve, gof] = fit(x,y,ft);
end

title(strcat(titleStr, sprintf("\n(%s)",fitfuncFullName(fitfunc))), ...
    "Interpreter","none");
if isnan(args.LegendColumns)
    numCol = ceil(numel(labels)/16);
else
    numCol = args.LegendColumns;
end
legend(labels, "Location", args.LegendPosition, "Interpreter","none", ...
    'NumColumns', numCol);
ylabel(ylabelStr);
xlabel(xlabelStr);
xlim([xMin xMax]);
ylim([yMin yMax]);
yticks([]);
hold off
set(gcf,'color','w');
box on

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end

