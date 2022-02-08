function FitResult = fit_T1(path, keywords, fitfunc, MarkerSize)
% Plot data files in given path based on given keywords. The keywords can 
% be as many as you want, and can be in any sequence.

arguments
    path string
    keywords string
    fitfunc string {mustBeMember(fitfunc,{'Exp1','Exp2','StrExp'})}= 'Exp1'
    MarkerSize double = 20
end

% Get the list of file names
fitfuncFullName = containers.Map( ...
    ["Exp1", "Exp2", "StrExp"], ...
    ["Mono-exponential", "Bi-exponential", "Stretched Exponential"]);

% Get the list of file names
keywords(end+1) = "T1PF";
keywords = unique(keywords);
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
labels = cell(2*length(files),1);

% Set the x and y axis label
xlabelStr = "Recovery time (\mus)";
ylabelStr = "Signal";

xMax = 0;
xMin = inf;
yMax = 0;
yMin = inf;
colorList = jet(length(files));
switch fitfunc
    case "Exp1"
        func = "c1 + c2*exp(-k*x)";
        num_k = 1; num_c = 2;
    case "Exp2"
        func = "c1 + c2*exp(-k1*x) + c3*exp(-k2*x)";
        num_k = 2; num_c = 3;
    case "StrExp"
        func = "c1 + c2*exp(-(k*x)^c3)";
        num_k = 1; num_c = 3;
end
k = zeros(length(files),num_k);
c = zeros(length(files),num_c);
rmse = zeros(size(files)); % Root mean squared error
FitResult = table(temperature, c, k, rmse);
for i = 1:length(files)
    f = files{i};
    [x, y] = eprload(fullfile(path, f));
    y = real(y);
    y = y/max(y); % normalize
    % change x axis
    DSCfile = strcat(f(1:end-3),"DSC");
    x = change_x_axis(numel(x), fullfile(path,DSCfile));
    x = x/1000;
    % fit
    % first use exponfit to get an estimate
    [k,c,yfit] = exponfit(x,y,1);
    switch fitfunc
        case "Exp1"
            fo = fitoptions(func,"StartPoint",[c k]);
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            FitResult.c(i,:) = [curve.c1 curve.c2];
            FitResult.k(i,:) = curve.k;
            FitResult.rmse(i) = gof.rmse; 
            yfit = curve(x);
        case "Exp2"
            bounds = [c/100 c(2)/100 k/100 k/100;...
                      c*100 c(2)*100 k*100 k*100];
            fo = fitoptions(func, ...
                "StartPoint",[c c(2) k k], ...
                "Lower",min(bounds,[],1), ...
                "Upper",max(bounds,[],1));
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            FitResult.c(i,:) = [curve.c1 curve.c2 curve.c3];
            FitResult.k(i,:) = [curve.k1 curve.k2];
            FitResult.rmse(i) = gof.rmse; 
            yfit = curve(x);
        case "StrExp"
            bounds = [c/100 0.1 k/100;...
                      c*100 10 k*100];
            fo = fitoptions(func, ...
                "StartPoint",[c 1 k], ...
                "Lower",min(bounds,[],1), ...
                "Upper",max(bounds,[],1));
            ft = fittype(func,"options",fo);
            [curve, gof] = fit(x,y,ft);
            FitResult.c(i,:) = [curve.c1 curve.c2 curve.c3];
            FitResult.k(i,:) = curve.k;
            FitResult.rmse(i) = gof.rmse; 
            yfit = curve(x);
    end
    scatter(x,y,MarkerSize, "MarkerEdgeColor",colorList(i,:));
    labels{2*i-1} = strcat(num2str(temperature(i))," K data");
    hold on
    plot(x,yfit, "Color", colorList(i,:));
    labels{2*i} = strcat(num2str(temperature(i))," K Fit");
    hold on
    xMin = min(xMin, min(x));
    xMax = max(xMax, max(x));
    yMin = min(yMin, min(real(y)));
    yMax = max(yMax, max(real(y)));
end
set(gca,'xscale','log');
title(strcat(titleStr, sprintf("\n(%s)",fitfuncFullName(fitfunc))), ...
    "Interpreter","none");
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


