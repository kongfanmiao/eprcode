function output = find_parameter(DSCfilePath, param)
% Find the parameters in DSC file
if ~endsWith(DSCfilePath, '.DSC')
    DSCfilePath = strcat(DSCfilePath, '.DSC');
end
fstr = fileread(DSCfilePath); % the whole file string
% Find the file string that matches the given pattern
paramStr = regexp(fstr, strcat(param, "\s+=\s+\d+\s+;"), "match");
paramStr = paramStr{:};
paramStr = split(paramStr);
output = str2double(paramStr{end-1});
end