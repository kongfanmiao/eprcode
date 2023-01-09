function save_T1(FitResult, path, fileName)
% Save T1 values. Always save in ns unit

temp = FitResult.Temperature;

try
    T1 = FitResult.T1; % mono or stretch exponential
catch
    try
        T1 = FitResult.T_long; % bi-exponential
    catch
        error('There is no T1 or T_long in the FitResult table');
    end
end

timeUnit = FitResult.Properties.VariableUnits{2};
switch timeUnit
    case 'ns'
        timeScaleFactor = 1;
    case '\mus'
        timeScaleFactor = 1e3;
    case 'ms'
        timeScaleFactor = 1e6;
    case 's'
        timeScaleFactor = 1e9;
end
T1 = T1*timeScaleFactor;
data = table(temp, T1);
data.Properties.VariableNames = {'Temperature', 'T1'};
data.Properties.VariableUnits = {'K', 'ns'};

filePath = fullfile(path, [fileName, '.csv']);
writecell(data.Properties.VariableNames, filePath);
writecell(data.Properties.VariableUnits, filePath, 'WriteMode','append');
writetable(data, filePath, 'WriteVariableNames',false, 'WriteMode','append');


end