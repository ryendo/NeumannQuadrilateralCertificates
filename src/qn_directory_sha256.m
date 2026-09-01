function digest = qn_directory_sha256(path)
% Hash a directory tree, excluding INTLAB's generated startup cache.

path=qn_canonical_path(path);
assert(isfolder(path),'Directory to hash does not exist: %s',path);
assert(~contains(path,''''), ...
    'Directory paths containing a single quote are unsupported.');
command=sprintf([ ...
    'cd ''%s'' && find . -type f ' ...
    '! -name ''Matlab_*_Intlab_Version_*.mat'' -print0 ' ...
    '| LC_ALL=C sort -z | xargs -0 sha256sum | sha256sum'],path);
[status,output]=system(command);
token=regexp(strtrim(output),'^([0-9a-f]{64})\s','tokens','once');
assert(status==0 && ~isempty(token),'Cannot hash directory %s.',path);
digest=token{1};
end
