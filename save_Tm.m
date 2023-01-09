function save_Tm(FitResult, path, fileName)
% Save Tm values. Always save in ns unit

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
Tm = Tm*timeScaleFactor;
data = table(temp, Tm);
data.Properties.VariableNames = {'Temperature', 'Tm'};
data.Properties.VariableUnits = {'K', 'ns'};

filePath = fullfile(path, [fileName, '.csv']);
writecell(data.Properties.VariableNames, filePath);
writecell(data.Properties.VariableUnits, filePath, 'WriteMode','append');
writetable(data, filePath, 'WriteVariableNames',false, 'WriteMode','append');


end