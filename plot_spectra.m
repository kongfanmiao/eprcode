function plot_spectra(varargin)
% Plot data files in given path based on given keywords. The keywords can 
% be as many as you want, and can be in any sequence.

% parse the input arguments
par = inputParser;

par.addRequired('path', @(x)isstring(x)|ischar(x));
par.addRequired('keywords', @isstring);
par.addParameter('Color', 'map');
par.addParameter('Offset', 0, @isnumeric);
par.addParameter('NormalizeField', false, @islogical);
par.addParameter('NormalizeSignal', false, @islogical);
par.addParameter('ScaleSignal', 1, @isnumeric);
par.addParameter('Shift', 0, @isnumeric);

par.KeepUnmatched = true;

parse(par, varargin{:});
path = par.Results.path;
keywords = par.Results.keywords;
color = par.Results.Color;
offset = par.Results.Offset;

% offset = 0.4;

% Get the list of file names
files = find_files(path, keywords(:));
% Extract the temperature string according to given pattern
tempStr = cellfun(@(s)regexp(s,"\d*\.*\d*K", "match"), files);
% Convert to double data type
temperature = cellfun(@(s)str2double(s(1:end-1)),tempStr);
% Sort the files according to temperature in ascending order
filesTable = table(files, temperature);
filesTable = sortrows(filesTable, "temperature");
temperature = filesTable.temperature;
files = filesTable.files;
% Convert cell array to matrix
keywordsList = horzcat(keywords(:)); % string array
% Create the figure title based on keywords
titleStr = join(keywordsList, ' ');
% Create figure legends
labels = cell(size(files));

% Set the x and y axis label
xlabelStr = "B (mT)";
ylabelStr = "Signal";

xMax = 0;
xMin = inf;
yMax = 0;
yMin = inf;
colorList = jet(length(files));

tmp = table();
for i = 1:length(files)
    f = files{i};
    [x, y, params] = eprload(fullfile(path, f));
    x = x/10; % use mT as unit
    y = real(y);
    
    % shorts per loop
    spp = params.ShotsPLoop;
    y = y/spp; % average all the shots

%     power = strip(params.Power);
%     power = power(1:end-2);
%     power = str2double(strip(power));
%     atten = strip(params.Attenuation);
%     atten = atten(1:end-2);
%     atten = str2double(strip(atten));
%     tmp{end+1,:} = [temperature(i), power, atten, max(y)];
   
    % Normalize the field
    if par.Results.NormalizeField
        % determine X or Q band
        mwfq = params.MWFQ;
        if (8e9 < mwfq) && (mwfq < 10e9) % X band, normalize to 9.5 GHz
            x = x*9.5e9/mwfq;
        elseif (33e9 < mwfq) && (mwfq < 36e9) % Q band, normalize to 34 GHz
            x = x*34e9/mwfq;
        end
    end          
    
    % Normalize the signal
    if par.Results.NormalizeSignal
        y = y/max(y);
    end
    
    % Offset of multiple lines
    y = y+ max(y) * offset;
    
    % Shift the line vertically
    y = y + max(y) * par.Results.Shift;

    % Scale the line verticalluy
    y = y*par.Results.ScaleSignal;

    if all(color == 'map')
        c = colorList(i,:);
    else
        c = color;
    end
    
    plot(x, y, "Color",c);
    labels{i} = strcat(num2str(temperature(i))," K");
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
ylim([yMin yMax]);
yticks([]);
hold off

fig = gcf;
if ~(fig.WindowStyle == "docked")
    set(fig,'position',[10,10,900,600]);
end

end


