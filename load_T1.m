function result = load_T1(path, keywords)
% load saved table

file = find_files(path,keywords,'FileType','csv', 'FullPath', true);
if numel(file)>1
    error('More than one file is found:')
    disp(file);
end
result = readtable(file{:});
result.Properties.VariableUnits = {'K', 'ns'};

end