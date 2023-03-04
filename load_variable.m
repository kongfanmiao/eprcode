function result = load_variable(path, keywords)
% load saved variable from .mat file

file = find_files(path,keywords, 'FileType','mat', 'FullPath', true);
if numel(file)>1
    error('More than one file is found:')
    disp(file);
end
result = load(file{:});

end