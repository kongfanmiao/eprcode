function plot_data(path, keywords)
% Plot data files in given path based on given keywords. The keywords can 
% be as many as you want, and can be in any sequence.

arguments
    path string
end
arguments (Repeating)
    keywords string
end

% Get the list of file names
files = find_files(path, keywords{:});
% Extract the temperature string according to given pattern
tempStr = cellfun(@(s)regexp(s,"\d+K", "match"), files);
% Convert to double data type
temperature = cellfun(@(s)str2double(s(1:end-1)),tempStr);
% Sort the files according to temperature in ascending order
filesTable = table(files, temperature);
filesTable = sortrows(filesTable, "temperature");
files = filesTable.files;
% Convert cell array to matrix
keywordsList = horzcat(keywords{:}); % string array
% Create the figure title based on keywords
titleStr = join(keywordsList, ' ');
% Create figure legends
labels = cellfun(@(s)s(5:end-4), files, "UniformOutput",false);

% Set the x and y axis label, as well as vertical offset scale when
% multiple lines are plot in same figure
if any(contains(keywordsList, 'EDFS'))
    xlabelStr = "B (mT)";
    ylabelStr = "Echo signal";
    offset = 0.4;
elseif any(contains(keywordsList, 'FIDFS'))
    xlabelStr = "B (mT)";
    ylabelStr = "FID signal";
    offset = 0.4;
elseif any(contains(keywordsList, 'T1PF'))
    xlabelStr = "Recovery time (\mus)";
    ylabelStr = "Echo signal";
    offset = 0;
elseif any(contains(keywordsList, 'Tm'))
    xlabelStr = strcat("2",char([0xD835 0xDF0F]), " (ns)");
    ylabelStr = "Echo signal";
    offset = 0;
end

xMax = 0;
xMin = inf;
yMax = 0;
yMin = inf;
colorList = jet(length(files));

for i = 1:length(files)
    f = files{i};
    [x, y] = eprload(fullfile(path, f));
    y = y/max(y); % normalize
    if any(keywordsList == 'T1PF')
        DSCfile = strcat(f(1:end-3),"DSC");
        x = change_x_axis(numel(x), fullfile(path,DSCfile));
        x = x/1e3; % change to um unit
        semilogx(x,y, "Color",colorList(i,:));
    else
        y = y+i*offset;
        plot(x, y, "Color",colorList(i,:));
    end
    hold on
    xMin = min(xMin, min(x));
    xMax = max(xMax, max(x));
    yMin = min(yMin, min(real(y)));
    yMax = max(yMax, max(real(y)));
end
title(titleStr, "Interpreter","none");
legend(labels, "Location","northeastoutside", "Interpreter","none");
ylabel(ylabelStr);
xlabel(xlabelStr);
xlim([xMin xMax]);
ylim([yMin yMax]);
yticks([]);
hold off

end


