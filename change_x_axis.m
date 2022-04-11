function xNew = change_x_axis(xOld, DSCfilePath)

% Change the x axis of T1 and Tm measurement
% xNew: one-dimensional data array
% DSCfilePath: path of the DSC file

arguments
    xOld double
    DSCfilePath char
end

N = numel(xOld);
if endsWith(DSCfilePath, '.DTA')
    DSCfilePath = strcat(DSCfilePath(1:end-3),'DSC');
end

p0 = find_parameter(DSCfilePath, "p0");
p2 = find_parameter(DSCfilePath, "p2");
d2 = find_parameter(DSCfilePath, "d2");
d3 = find_parameter(DSCfilePath, "d3");
d30 = find_parameter(DSCfilePath, "d30");
% the i-th value of xNew = d3+d30*(1+2+...+i-1)+0.5*p9-p2
xNew = d3+d30*(1:N)'.*(0:N-1)'/2+0.5*p0-p2;
end

function output = find_parameter(DSCfilePath, param)
% Find the parameters in DSC file
fstr = fileread(DSCfilePath); % the whole file string
% Find the file string that matches the given pattern
paramStr = regexp(fstr, strcat(param, "\s+=\s+\d+\s+;"), "match");
paramStr = paramStr{:};
paramStr = split(paramStr);
output = str2double(paramStr{end-1});
end