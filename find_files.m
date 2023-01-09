% find_files    Find files in given path based on give keyword
%   find_files(path, {'keyword1', 'keyword2'})
%   find_files(path, {'keyword1', 'keyword2'}, 'FileType', 'DTA')
%   find_files(path, {'keyword1', 'keyword2'}, 'FileType', 'csv')
%   find_files(path, {'keyword1', 'keyword2'}, 'FullPathOutput', true)
%
%   Input:
%       path        directory to find files
%       keywords    keywords based on which to find files
%                   - string: "80K" (if there's only one keyword)
%                   - character vector: '80K' (if there's only one keyword)
%                   - string array: ["80K", "3480G"]
%                   - cell array of character vectors: {'80K', '3480G'}
%       FileType    optional, can be any type, like .DTA, .csv
%       FullPathOutput
%                   optional, if true spit out the absolute path of the
%                   file(s). If false spit out only the file names
%
%   Output:
%       cell array containing the files


function outputFiles = find_files(path, keywords, varargin)

par = inputParser;
par.addParameter('FileType', "DTA", @(x)isstring(x)|ischar(x));
par.addParameter('FullPathOutput', false, @islogical);
par.KeepUnmatched = true;

parse(par, varargin{:});
args = par.Results;

% Add file fype to keywords
keywords = string(keywords);
keywords(end+1) = strcat(".", args.FileType); 

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

if args.FullPathOutput
    outputFiles = cellfun(@(file)fullfile(path,file), outputFiles, ...
        'UniformOutput',false);
end

% call error if no files are found
if isempty(outputFiles)
    error("I couldn't find any file that contains all of the keywords")
end

end