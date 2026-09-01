function [commit,dirty] = qn_git_state(path,ignored_untracked_prefixes)
% Return the checked-out Git commit and whether tracked files are modified.

if nargin<2, ignored_untracked_prefixes={}; end

assert(ischar(path) || (isstring(path) && isscalar(path)), ...
    'Git path must be a character vector or scalar string.');
path=char(path);
assert(isfolder(path),'Git path does not exist: %s',path);
assert(~contains(path,''''),'Git paths containing a single quote are unsupported.');
quoted=['''' path ''''];
[status,commit]=system(sprintf('git -C %s rev-parse HEAD 2>/dev/null',quoted));
commit=strtrim(commit);
assert(status==0 && ~isempty(regexp(commit,'^[0-9a-f]{40}$','once')), ...
    'Cannot determine a full Git commit for %s.',path);
[status,output]=system(sprintf( ...
    'git -C %s status --porcelain --untracked-files=all 2>/dev/null',quoted));
assert(status==0,'Cannot determine Git status for %s.',path);
dirty=false;
lines=regexp(strtrim(output),'\r?\n','split');
for k=1:numel(lines)
    line=char(lines{k});
    if isempty(line), continue; end
    if startsWith(line,'?? ')
        relative_path=line(4:end);
        ignored=any(cellfun(@(prefix)startsWith(relative_path,prefix), ...
            ignored_untracked_prefixes));
        if ignored, continue; end
    end
    dirty=true;
    break
end
end
