function digest = qn_sha256_file(path)
% Compute a file SHA-256 digest with the platform sha256sum utility.

assert(ischar(path) || (isstring(path) && isscalar(path)), ...
    'Hash path must be a character vector or scalar string.');
path=char(path);
assert(isfile(path),'File to hash does not exist: %s',path);
assert(~contains(path,''''),'Hash paths containing a single quote are unsupported.');
[status,output]=system(sprintf('sha256sum -- ''%s'' 2>/dev/null',path));
token=regexp(strtrim(output),'^([0-9a-f]{64})\s','tokens','once');
assert(status==0 && ~isempty(token),'Cannot compute SHA-256 for %s.',path);
digest=token{1};
end
