function canonical = qn_canonical_path(path)
% Resolve a filesystem path without requiring a particular spelling.

assert(ischar(path) || (isstring(path) && isscalar(path)), ...
    'Path must be a character vector or scalar string.');
file=javaObject('java.io.File',char(path));
canonical=char(file.getCanonicalPath());
end
