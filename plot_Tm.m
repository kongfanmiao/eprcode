function FitResult = plot_Tm(path, keywords, fit, nExps)
% Plot data files in given path based on given keywords. The keywords can 
% be as many as you want, and can be in any sequence.

arguments
    path string
    keywords string
    fit = false
    nExps = 1
end

% Get the list of file names
keywords(end+1) = "Tm";
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
titleStr = join(keywordsList, " ");
% Create figure legends
labels = {};

% Set the x and y axis label
xlabelStr = strcat("2",char([0xD835 0xDF0F]), " (ns)");
ylabelStr = "Signal";

xMax = 0;
xMin = inf;
yMax = 0;
yMin = inf;
colorList = jet(length(files));

k = zeros(length(files),nExps);
c = zeros(length(files),nExps+1);
FitResult = table(temperature, k, c);
for i = 1:length(files)
    f = files{i};
    [x, y] = eprload(fullfile(path, f));
    y = real(y);
    y = y/max(y); % normalize
    % fit
    if fit
        [k,c,yfit] = exponfit(x,y,nExps);
        FitResult.k(i,:) = k;
        FitResult.c(i,:) = c;
        scatter(x,y,20, "MarkerEdgeColor",colorList(i,:));
        labels{end+1} = strcat(num2str(temperature(i))," K");
        hold on
        plot(x,yfit, "Color", colorList(i,:));
        labels{end+1} = strcat(num2str(temperature(i))," K Fit");
    else
        plot(x,y, "Color",colorList(i,:));
        labels{end+1} = strcat(num2str(temperature(i))," K");
    end
    hold on
    xMin = min(xMin, min(x));
    xMax = max(xMax, max(x));
    yMin = min(yMin, min(real(y)));
    yMax = max(yMax, max(real(y)));
end
% set(gca,'xscale','log');
title(titleStr, "Interpreter","none");
legend(labels, "Location","northeastoutside", "Interpreter","none");
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


