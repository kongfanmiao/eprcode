function result = load_csv(path, keywords)
% load saved table from csv file. It can be T1, Tm, or simulated spectrum.
% The second row is the unit

file = find_files(path,keywords,'FileType','csv', 'FullPath', true);
if numel(file)>1
    error('More than one file is found:')
    disp(file);
end
% read from the third row
file = file{:};
opts = detectImportOptions(file);
opts.VariableNamesLine = 1; % row number which has variable names
opts.DataLine = 3; % row number from which the actual data starts
opts.VariableUnitsLine = 2; % row number which has units specified
result = readtable(file, opts);

end