function outputFiles = find_files(path, keywords)
% Find files in given path based on give keywords. The keywords can be as 
% many as you want, and can be in any sequence. They keywords list should be
% a string array.
arguments
    path string
    keywords string % string array
end

% check what file type to read
if ~any(contains(keywords, "DTA") | contains(keywords, "DSC"))
    keywords(end+1) = "DTA"; % by default read DTA file, if no file 
                                  % type is provided
end

% if the temperature keywords is provides, prepend a '_' before searching
% the file
idxK = contains(keywords, "K");
if any(idxK) % if constains a temperature string
    keywords(idxK) = strcat("_",keywords(idxK));
end

% Now search the files that contain all of the keywords
allFiles = {dir(path).name}';
indices = cellfun(@(s)contains(allFiles,s), keywords, ...
    "UniformOutput",false);
finIdx = all(horzcat(indices{:}),2);
outputFiles = allFiles(finIdx); % one column cell array

% call error if no files are found
if isempty(outputFiles)
    error("I couldn't find any file that contains all of the keywords")
end

end