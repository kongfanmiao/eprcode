function save_T1(FitResult, path, fileName)

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

data = table(temp, T1);
data.Properties.VariableNames = {'Temperature', 'T1'};
data.Properties.VariableUnits = {'K', timeUnit};

filePath = fullfile(path, [fileName, '.csv']);
writecell(data.Properties.VariableNames, filePath);
writecell(data.Properties.VariableUnits, filePath, 'WriteMode','append');
writetable(data, filePath, 'WriteVariableNames',false, 'WriteMode','append');


end