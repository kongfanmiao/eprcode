function save_simulation(B, spc, fitResult, path, dataInfo)

fitSpc = reshape(fitResult.fit, [], 1);
fitSys = fitResult.argsfit{1};
data = table(B, spc, fitSpc);
data.Properties.VariableNames = {'B', 'Signal', 'Simulated signal'};
data.Properties.VariableUnits = {'mT', 'arb. un.', 'arb. un.'};

spcFileName = [char(dataInfo) '_SimulatedSpectrum'];
spcFilePath = fullfile(path, [spcFileName '.csv']);
writecell(data.Properties.VariableNames, spcFilePath);
writecell(data.Properties.VariableUnits, spcFilePath, 'WriteMode','append');
writetable(data, spcFilePath, 'WriteVariableNames',false, 'WriteMode','append');

sysFileName = [char(dataInfo) '_SimulatedSystem'];
sysFilePath = fullfile(path, [sysFileName '.xml']);
if numel(fitSys) == 1
    writestruct(fitSys, sysFilePath);
else
    for i = 1:numel(fitSys)
        fitSysNew.(sprintf("Sys%d",i)) = fitSys{i};
    end
    writestruct(fitSysNew, sysFilePath);
end

end