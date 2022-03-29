function [x,y,params] = load_file(path, keywords)
% load data base on give keywords
arguments
    path string
    keywords string % string array
end

file = find_files(path, keywords);
if length(file) ~= 1
    disp(file);
    error("More than one file found. Please revise your keywords");
end
file = file{:};
[x,y,params] = eprload(fullfile(path,file));
end