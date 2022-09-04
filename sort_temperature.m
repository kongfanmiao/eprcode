function filesTable = sort_temperature(path, keywords)

% Find all the files in a give path based on keywords
% Then sort the files by temperature in ascending order. There should be
% temperature string such as **K in file name

Files = find_files(path, keywords);
% Extract the temperature string according to given pattern
tempStr = cellfun(@(s)regexp(s,"\d*\.*\d*K", "match"), Files);
% Convert to double data type
Temperature = cellfun(@(s)str2double(s(1:end-1)),tempStr);
% Sort the files according to temperature in ascending order
filesTable = table(Files, Temperature);
filesTable = sortrows(filesTable, "Temperature");

end