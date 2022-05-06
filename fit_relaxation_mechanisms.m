% fit_relaxation_mechanisms     fit 1/T1 ~ T
%
%   fit_relaxation_mechanisms(FitResult, {'direct', 'raman'})
%   fit_relaxation_mechanisms(T, T1, {'direct', 'raman'})
%
%   Input
%       FitResult   table containing T1 fitting result
%       T           temperature array
%       T1          T1 array
%       (either provide FitResult or T and T1)
%       Mechanims   relaxation mechanism to use
%                   {'direct', 'raman'}
%       MwFreq      provide microwave frequency if using thermal process,
%                   by default in GHz unit. (optional)
%       LogScale    plot x and y axis in log scale
%       PlotIndividualProcess
%                   plot each summand of relaxation processes
%       TimeUnit    time unit of the T1 values. This is required if you
%                   provide T and T1 as array but not table
%       StartingPoint
%                   starting point of the fitting
%       LowerBound  lower bound for the fitting
%       UpperBound  upper bound for the fitting
%
%   Output
%       curve       fitted curve
%       gof         goodness of fitting
%


function varargout = fit_relaxation_mechanisms(varargin)

allMechanisms = {'direct', 'raman', 'orbach', 'local', 'thermal', 'power'};


par = inputParser;
par.KeepUnmatched = true;
par.addParameter('MwFreq', nan, @isnumeric);
par.addParameter('LogScale', false, @islogical);
par.addParameter('PlotIndividualProcess', false, @islogical);
par.addParameter('TimeUnit', nan, @(x)isstring(x)|ischar(x));
par.addParameter('StartingPoint', nan, @isnumeric);
par.addParameter('LowerBound', nan, @isnumeric);
par.addParameter('UpperBound', nan, @isnumeric);


% PARSE ARGUMENTS
if isa(varargin{1}, 'table')
    FitResult = varargin{1};
    T = FitResult.Temperature;
    try
        T1 = FitResult.T1; % mono or stretch exponential
    catch
        try
            T1 = FitResult.T_long; % bi-exponential
        catch
            error('There is no T1 or T_long in the FitResult table');
        end
    end
    T1inv = 1./T1;
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
        T1inv = 1./T1;
        mechanisms = string(varargin{3});
        mechanisms = arrayfun(@(s)validatestring(s, allMechanisms), ...
            mechanisms, 'UniformOutput', false);
        par.parse(varargin{4:end});
        unit = par.Results.timeUnit;
    catch
        error(['If the first argument is numeric, the second argument ' ...
            'must be also numeric and with the same length']);
    end
end
logT1inv = log10(T1inv);

args = par.Results;
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

% DEFINE FITTING FUNCTIONS FOR DIFFERENT RELAXATION PROCESSES
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
power = @(A_pow, n, T) A_pow*T.^n;

% CONSTRUCT FINAL FITTING FUNCTION
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

% GUESS STARTING POINT AND DETERMINE BOUNDS
% Direct process
A_dir = mean(T1inv)/mean(T);
A_dirMin = A_dir/1e4; A_dirMax = A_dir*1e4;
% Raman process
A_ram = A_dir;
A_ramMin = A_dir/1e4; A_ramMax = A_dir*1e4;
theta_D = 100; % K
theta_DMin = 0; theta_DMax = 1e5;
% Local process
A_loc = A_dir;
A_locMin = A_dir/1e4; A_locMax = A_dir*1e4;
Delta_loc = 100; % K
Delta_locMin = 0; Delta_locMax = 1e5;
% Orbach process
Delta_orb = 100; % K
Delta_orbMin = 0; Delta_orbMax = 1e4;
A_orb = mean(T1inv)/Delta_orb^3;
A_orbMin = A_orb/1e4; A_orbMax = A_orb*1e4;
% Thermal process
if contains('thermal', mechanisms)
    A_therm = omega^3;
    A_thermMin = A_therm/1e4; A_thermMax = A_therm*1e4;
    tau_0 = 1;
    tau_0Min = tau_0/1e4; tau_0Max = tau_0*1e4;
    E_a = 10;
    E_aMin = 0; E_aMax = 1e4;
end
% Power
A_pow = A_dir;
A_powMin = A_pow/1e4; A_powMax = A_pow*1e4; 
n = 2;
nMin = 1; nMax = 7;

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

% FIT
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

% PLOT
scatter(T, T1inv, 'o', 'blue');
hold on
plot(T, 10.^curve(T), '-k');
legend('',strjoin(mechanisms,'+'));

if args.PlotIndividualProcess
    for i = 1:numel(mechanisms)
        func = eval(mechanisms{i});
        coefsName = coefsGroup{i};
        coefsValue = cellfun(@(x){curve.(x)}, coefsName);
        plot(T, func(coefsValue{:}, T), '--');
    end
    legend('',strjoin(mechanisms,'+'),mechanisms{:}, "location", 'best');
end
hold off

xlabel('Temperature (K)');
ylabel(sprintf('{T_1}^{-1} (%s^{-1})', unit));
title({'Fit spin-lattice relaxation mechanisms', ...
       ['(', strjoin(mechanisms, ', '), ')']});

if args.LogScale
    set(gca, 'XScale', 'log', 'YScale', 'log');
end

set(gcf,'color','w');
box on

fig = gcf;
if fig.WindowStyle ~= "docked"
    set(fig,'position',[10,10,900,600]);
end

end
