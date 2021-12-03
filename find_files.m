function outputFiles = find_files(path, keywords)
% Find files in given path based on give keywords. The keywords can be as 
% many as you want, and can be in any sequence
arguments
    path string
end
arguments (Repeating)
    keywords string
end

keywordsList = horzcat(keywords{:}); % string array
% check what file type to read
if ~any(contains(keywordsList, "DTA") | contains(keywordsList, "DSC"))
    keywordsList(end+1) = "DTA"; % by default read DTA file, if no file 
                                  % type is provided
end

% if the temperature keywords is provides, prepend a '_' before searching
% the file
idxK = contains(keywordsList, "K");
if any(idxK) % if constains a temperature string
    keywordsList(idxK) = strcat("_",keywordsList(idxK));
end

% Now search the files that contain all of the keywords
allFiles = {dir(path).name}';
indices = cellfun(@(s)contains(allFiles,s), keywordsList, ...
    "UniformOutput",false);
finIdx = all(horzcat(indices{:}),2);
outputFiles = allFiles(finIdx); % one column cell array

% call error if no files are found
if isempty(outputFiles)
    error("I couldn't find any file that contains all of the keywords")
end

end