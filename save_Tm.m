function save_Tm(FitResult, path, fileName)

temp = FitResult.Temperature;

try
    Tm = FitResult.Tm; % mono or stretch exponential
catch
    try
        Tm = FitResult.T_long; % bi-exponential
    catch
        error('There is no Tm or T_long in the FitResult table');
    end
end

timeUnit = FitResult.Properties.VariableUnits{2};

data = table(temp, Tm);
data.Properties.VariableNames = {'Temperature', 'Tm'};
data.Properties.VariableUnits = {'K', timeUnit};

filePath = fullfile(path, [fileName, '.csv']);
writecell(data.Properties.VariableNames, filePath);
writecell(data.Properties.VariableUnits, filePath, 'WriteMode','append');
writetable(data, filePath, 'WriteVariableNames',false, 'WriteMode','append');


end