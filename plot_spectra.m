% plot_spectra    Plot the field sweep spectra
%
%   plot_spectra(path, keywords)
%   plot_spectra(path, keywords, 'Color','k')
%   plot_spectra(path, keywords, 'Offset', 0)
%   plot_spectra(path, keywords, 'NormalizeField', true)
%   plot_spectra(path, keywords, 'NormalizeSignal', true)
%   plot_spectra(path, keywords, 'ScaleSignal', 2)
%   plot_spectra(path, keywords, 'Shift', 0.5)
%   plot_spectra(path, keywords, 'FieldUnit', 'G')
%
%   Input:
%       path        directory to find files
%       keywords    keywords based on which to find files
%                   - string: "EDFS" (if there's only one keyword)
%                   - character vector: 'EDFS' (if there's only one keyword)
%                   - string array: ["80K", "EDFS"]
%                   - cell array of character vectors: {'80K', 'EDFS'}
%       Color       color of the spectra, use jet colormap by default, you
%                   can specify a color and it applies to all lines
%       Offset      apply a vertical offset if there are multiple lines. By
%                   default it's zero
%       NormalizeField    
%                   normalize the frequency and scale the field
%                   accordingl. X band is normalilzed to 9.5 GHz, Q band to
%                   34 GHz. It's false by default.
%       NormailzeSignal
%                   normalize the signal magnitude to 1 for all the lines.
%                   It's false by default
%       ScaleSignal stretch or compress the signal magnitude vertically.
%                   The default value is 1
%       Shift       shift the signal vertically. Default is 0.
%       FieldUnit   unit of field, 'mT' by default
%
%   Example:
%       plot_spectra('./Data',["EDFS","5K"] 'color', 'k')
%

function plot_spectra(path, keywords, varargin)

% parse the input arguments
par = inputParser;

par.addParameter('Color', 'k');
par.addParameter('Offset', 0, @isnumeric);
par.addParameter('NormalizeField', false, @islogical);
par.addParameter('NormalizeSignal', false, @islogical);
par.addParameter('ScaleSignal', 1, @isnumeric);
par.addParameter('Shift', 0, @isnumeric);
par.addParameter('FieldUnit', 'mT', @(x)any(validatestring(x, {'mT', 'G'})));

par.KeepUnmatched = true;

parse(par, varargin{:});
args = par.Results;

color = args.Color;
offset = args.Offset;

keywords = string(keywords); 
filesTable = sort_temperature(path, keywords);
Temperature = filesTable.Temperature;
Files = filesTable.Files;
numFiles = length(Files);

% Create the figure title based on keywords
titleStr = join(keywords, ' ');
% Create figure legends
labels = cell(size(Files));

% Set the x and y axis label
xlabelStr = sprintf("B (%s)", args.FieldUnit);
ylabelStr = "Signal (arb. u.)";

xMax = 0;
xMin = inf;
yMax = 0;
yMin = inf;
colorList = jet(numFiles);

for i = 1:numFiles
    f = Files{i};
    [x, y, params] = eprload(fullfile(path, f));
    if isequal(args.FieldUnit, 'mT')
        x = x/10;
    end
    y = real(y);
    % average all the shots
    spp = params.ShotsPLoop; % shorts per loop
    y = y/spp;
%     power = strip(params.Power);
%     power = power(1:end-2);
%     power = str2double(strip(power));
%     atten = strip(params.Attenuation);
%     atten = atten(1:end-2);
%     atten = str2double(strip(atten));
%     tmp{end+1,:} = [Temperature(i), power, atten, max(y)];
   
    % Normalize the field
    if args.NormalizeField
        % determine X or Q band
        mwfq = params.MWFQ;
        if (8e9 < mwfq) && (mwfq < 10e9) % X band, normalize to 9.5 GHz
            x = x*9.5e9/mwfq;
        elseif (33e9 < mwfq) && (mwfq < 36e9) % Q band, normalize to 34 GHz
            x = x*34e9/mwfq;
        end
    end          
    
    % Normalize the signal
    if args.NormalizeSignal
        y = y/max(y);
    end
    
    % Offset of multiple lines
    y = y+ max(y) * offset;
    
    % Shift the line vertically
    y = y + max(y) * args.Shift;

    % Scale the line verticalluy
    y = y*args.ScaleSignal;

    if isequal(color,'map')
        c = colorList(i,:);
    else
        c = color;
    end
    
    plot(x, y, "Color",c);
    labels{i} = strcat(num2str(Temperature(i))," K");
    hold on
    xMin = min(xMin, min(x));
    xMax = max(xMax, max(x));
    yMin = min(yMin, min(real(y)));
    yMax = max(yMax, max(real(y)));
end

% tmp.Properties.VariableNames = {'Temperature', 'Power', 'Attenuation', 'Max'};
% disp(tmp)

title(titleStr, "Interpreter","none");
legend(labels, "Location","best", "Interpreter","none");
ylabel(ylabelStr);
xlabel(xlabelStr);
xlim([xMin xMax]);
ylim([yMin yMax*1.05]);
yticks([]);
hold off
set(gcf,'color','w');
box on

% fig = gcf;
% if fig.WindowStyle ~= "docked"
%     set(fig,'position',[10,10,900,600]);
% end

end

