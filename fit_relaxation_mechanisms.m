% fit_relaxation_mechanisms     fit 1/T1 ~ T
%
%   fit_relaxation_mechanisms(FitResult, {'direct', 'raman'})
%   fit_relaxation_mechanisms(T, T1, {'direct', 'raman'})
%
%   Input
%       FitResult   table containing T1 fitting result
%       T           temperature array
%       T1          T1 array, will be converted to seconds
%       (either provide FitResult or T and T1)
%       Mechanism   relaxation mechanism to use
%                   {'direct', 'raman'}
%       MwFreq      provide microwave frequency if using thermal process,
%                   by default in Hz unit. (optional)
%       LogScale    plot x and y axis in log scale
%       PlotIndividualProcess
%                   plot each summand of relaxation processes
%       TimeUnit    time unit of the T1 values. This is required if you
%                   provide T and T1 as array but not table
%       StartingPoint
%                   starting point of the fitting
%       LowerBound  lower bound for the fitting
%       UpperBound  upper bound for the fitting
%       SearchRange search range of parameters
%       MarkerSize  marker size
%       filled      filled markers
%
%   Output
%       curve       fitted curve
%       gof         goodness of fitting
%


function varargout = fit_relaxation_mechanisms(varargin)

allMechanisms = {'direct', 'raman', 'orbach', 'local', 'thermal', ...
    'cosech', 'power', 'cross'};


par = inputParser;
par.KeepUnmatched = true;
par.addParameter('MwFreq', nan, @isnumeric);
par.addParameter('LogScale', false, @islogical);
par.addParameter('PlotIndividualProcess', false, @islogical);
par.addParameter('TimeUnit', nan, @(x)isstring(x)|ischar(x));
par.addParameter('StartingPoint', nan, @isnumeric);
par.addParameter('LowerBound', nan, @isnumeric);
par.addParameter('UpperBound', nan, @isnumeric);
par.addParameter('SearchRange', 1e3, @isnumeric);
par.addParameter('MarkerSize', 100, @isnumeric);
par.addParameter('filled', true, @islogical);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse arguments
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(varargin{1}, 'table')
    FitResult = varargin{1};
    T = FitResult.Temperature;
    T1 = nan;
    if ismember('T1', FitResult.Properties.VariableNames)
        T1 = FitResult.T1; % mono or stretch exponential
    elseif ismember('T_long', FitResult.Properties.VariableNames)
        T1 = FitResult.T_long; % bi-exponential
    else
        error('There is no T1 or T_long in the FitResult table');
    end
    unit = FitResult.Properties.VariableUnits{2};
    mechanisms = string(varargin{2});
    mechanisms = arrayfun(@(s)validatestring(s, allMechanisms), ...
        mechanisms, 'UniformOutput', false);
    par.parse(varargin{3:end});
elseif isa(varargin{1}, 'numeric')
    try
        assert(isa(varargin{2}, 'numeric'));
        assert(numel(varargin{1})==numel(varargin{2}));
        T = varargin{1};
        T1 = varargin{2};
        mechanisms = string(varargin{3});
        mechanisms = arrayfun(@(s)validatestring(s, allMechanisms), ...
            mechanisms, 'UniformOutput', false);
        par.parse(varargin{4:end});
        unit = par.Results.TimeUnit;
    catch
        error(['If the first argument is numeric, the second argument ' ...
            'must be also numeric and with the same length']);
    end
end
% convert T1 to seconds. So that 1/T1 will be a large value
switch unit
    case 'ns'
        timeScaleFactor = 1e9;
    case '\mus'
        timeScaleFactor = 1e6;
    case 'ms'
        timeScaleFactor = 1e3;
    case 's'
        timeScaleFactor = 1;
end
T1 = T1/timeScaleFactor;
T1inv = 1./T1;
logT1inv = log10(T1inv);

args = par.Results;
SearchRange = args.SearchRange;
% Ask for microwave frequency if use thermal process
if contains('thermal', mechanisms)
    if ~isnan(args.MwFreq)
        omega = args.MwFreq;
    else
        error(['You should provide the microwave frequency if using ' ...
            'thermal process']);
    end
    % remove the 1e9 of frequenc, use ns as unit for correlation time tau
    omega = omega*2*pi/1e9; % rad/s
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define fitting functions for different relaxation processes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit log(1/T1) ~ log(T)
direct = @(A_dir, T) A_dir*T;
raman = @(A_ram, theta_D, T) ...
    arrayfun(@(xx)A_ram*xx.^9*integral( ...
    @(x)(x.^8.*exp(x))./((exp(x)-1).^2), 0, 1/xx), T/theta_D);
% Delta_loc is in unit of K
local = @(A_loc, Delta_loc, T) ...
    A_loc*(exp(Delta_loc./T))./(exp(Delta_loc./T)-1).^2;
% Delta_orb is in unit of K
orbach = @(A_orb, Delta_orb, T) (A_orb*Delta_orb^3)./(exp(Delta_orb./T)-1);
thermal = @(A_therm, tau_0, E_a, T) ...
    A_therm*(2*tau_0*exp(E_a./T))./(1+omega^2*(tau_0*exp(E_a./T)).^2);
cosech = @(A_tls, Delta_tls, T) ...
    A_tls*2./(exp(Delta_tls./T)-exp(-Delta_tls./T));
power = @(A_pow, n, T) A_pow*T.^n;
cross = @(A_cro, Delta_cro, T) A_cro./(1+exp(Delta_cro./T));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct final fitting function 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coefs = {};
coefsGroup = cell(1,numel(mechanisms));
funcStr = {};
for i = 1:numel(mechanisms)
    j1 = 0;
    j2 = 0;
    mechStr = mechanisms{i};
    fullFuncStr = functions(eval(mechStr)).function; % character vector
    for j = 1:numel(fullFuncStr)
        c = fullFuncStr(j);
        if isequal(c, '(')
            j1 = j;
        end
        if isequal(c, ')')
            j2 = j;
            break
        end
    end
    tmpCoef = split(fullFuncStr(j1+1:j2-3), ',')';
    coefsGroup{i} = tmpCoef;
    coefs = [coefs tmpCoef];
    funcStr = [funcStr fullFuncStr(j2+1:end)];
end
coefStr = ['@(', strjoin(coefs,','), ',T)'];
funcStr = ['log10(', strjoin(funcStr, '+') ,')'];
fitfuncStr = [coefStr, funcStr];
fitfunc = eval(fitfuncStr);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guess starting point and determine bounds
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct process
A_dir = mean(T1inv)/mean(T);
A_dirMin = A_dir/SearchRange; A_dirMax = A_dir*SearchRange;
% Raman process
A_ram = A_dir;
A_ramMin = A_dir/SearchRange; A_ramMax = A_dir*SearchRange;
theta_D = 100; % K
theta_DMin = 10; theta_DMax = 1e4;
% Local process
A_loc = A_dir;
A_locMin = A_dir/SearchRange; A_locMax = A_dir*SearchRange;
Delta_loc = 100; % K
Delta_locMin = 10; Delta_locMax = 1e5;
% Orbach process
Delta_orb = 100; % K
Delta_orbMin = 10; Delta_orbMax = Delta_orb*SearchRange;
A_orb = mean(T1inv)/Delta_orb^3;
A_orbMin = A_orb/SearchRange; A_orbMax = A_orb*SearchRange;
% Thermal process
if contains('thermal', mechanisms)
    A_therm = omega^3;
    A_thermMin = A_therm/SearchRange; A_thermMax = A_therm*SearchRange;
    tau_0 = 1;
    tau_0Min = tau_0/SearchRange; tau_0Max = tau_0*SearchRange;
    E_a = 10;
    E_aMin = 0; E_aMax = SearchRange;
end
% cosech
A_tls = mean(T1inv); Delta_tls = 50;
A_tlsMin = A_tls/SearchRange; A_tlsMax = A_tls*SearchRange;
Delta_tlsMin = 1; Delta_tlsMax = Delta_tls*SearchRange;
% Power
A_pow = A_dir;
A_powMin = A_pow/SearchRange; A_powMax = A_pow*SearchRange; 
n = 2;
nMin = 1; nMax = 7;
% cross
A_cro = mean(T1inv); Delta_cro = 10;
A_croMin = A_cro/SearchRange; A_croMax = A_cro*SearchRange;
Delta_croMin = 1; Delta_croMax = Delta_cro*SearchRange;

% Define starting point and upper, lower bound
startingPoint = zeros(1,numel(coefs));
bounds = zeros(2, numel(coefs));
for i = 1:numel(coefs)
    coef = coefs{i};
    startingPoint(i) = eval(coef);
    bounds(1,i) = eval([coef, 'Min']);
    bounds(2,i) = eval([coef, 'Max']);
end
lowerBound = min(bounds,[],1);
upperBound = max(bounds,[],1);

if ~isnan(args.StartingPoint)
    startingPoint = args.StartingPoint;
end
if ~isnan(args.LowerBound)
    lowerBound = args.LowerBound;
end
if ~isnan(args.UpperBound)
    upperBound = args.UpperBound;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ft = fittype(fitfunc, 'independent', 'T');
[curve, gof] = fit(T, logT1inv, ft, ...
    'StartPoint', startingPoint, ...
    'Upper', upperBound, ...
    'Lower', lowerBound);
if nargout == 0
    varargout = {};
elseif nargout == 1
    varargout{:} = curve;
elseif nargout == 2
    varargout{:} = [curve, gof];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

additionalArgs = {};
if args.filled
    additionalArgs = {'filled'};
end
% Plot data in seconds
scatter(T, T1inv, args.MarkerSize, 'o', 'blue', additionalArgs{:});
hold on
Tnew = min(T)*0.98:max(T)*1.02;
plot(Tnew, 10.^curve(Tnew), '-k', 'LineWidth', 1.5);
legend('',strjoin(mechanisms,'+'), 'location', 'best');

if args.PlotIndividualProcess
    for i = 1:numel(mechanisms)
        func = eval(mechanisms{i});
        coefsName = coefsGroup{i};
        coefsValue = cellfun(@(x){curve.(x)}, coefsName);
        plot(Tnew, func(coefsValue{:}, Tnew), '--');
    end
    legend('',strjoin(mechanisms,'+'),mechanisms{:}, "location", 'best');
end
hold off

xlabel('Temperature (K)');
ylabel('{T_1}^{-1} (s^{-1})');
ylim([0 max(T1inv)*1.1]);
if args.LogScale
    ylim([min(T1inv)/2 max(T1inv)*2]);
end
title({'Fit spin-lattice relaxation mechanisms', ...
       ['(', strjoin(mechanisms, ', '), ')']});

if args.LogScale
    set(gca, 'XScale', 'log', 'YScale', 'log');
end

set(gcf,'color','w');
box on

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end
