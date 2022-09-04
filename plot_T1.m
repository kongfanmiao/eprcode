function plot_T1(path, keywords)
% Plot data files in given path based on given keywords. The keywords can 
% be as many as you want, and can be in any sequence.

arguments
    path char
    keywords string
end

% Get the list of file names
% keywords = string(keywords);
keywords(end+1) = "T1PF";
keywords = unique(keywords);
filesTable = sort_temperature(path, keywords);
Temperature = filesTable.Temperature;
Files = filesTable.Files;
numFiles = length(Files);

% Create the figure title based on keywords
titleStr = join(keywords, " ");
% Create figure legends
labels = cell(size(Files));

% Set the x and y axis label
xlabelStr = "Recovery time (\mus)";
ylabelStr = "Signal (arb. u.)";

xMax = 0;
xMin = inf;
yMax = 0;
yMin = inf;
colorList = jet(numFiles);

for i = 1:numFiles
    f = Files{i};
    [x, y] = eprload(fullfile(path, f));
    y = real(y);
    y = y/max(y); % normalize
    % change x axis
    x = change_x_axis(x, fullfile(path,f));
    x = x/1e3; % change to us unit
    plot(x,y, "Color",colorList(i,:));
    labels{i} = strcat(num2str(Temperature(i))," K");
    hold on
    xMin = min(xMin, min(x));
    xMax = max(xMax, max(x));
    yMin = min(yMin, min(real(y)));
    yMax = max(yMax, max(real(y)));
end
set(gca,'xscale','log');
title(titleStr, "Interpreter","none");
legend(labels, "Location","best", "Interpreter","none");
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


