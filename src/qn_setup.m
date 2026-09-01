function repo_root=qn_setup(intlab_root,veigs_root)
% Add the proof code and its two verified numerical dependencies.

repo_root=fileparts(fileparts(mfilename('fullpath')));
if nargin>=1 && ~isempty(intlab_root), addpath(intlab_root); end
if nargin>=2 && ~isempty(veigs_root), addpath(veigs_root); end
addpath(fullfile(repo_root,'src'),fullfile(repo_root,'tests'));
assert(exist('startintlab','file')==2,'INTLAB startup file is unavailable.');
intlab_start_file=which('startintlab');
if nargin>=1 && ~isempty(intlab_root)
    assert(strcmp(qn_canonical_path(fileparts(intlab_start_file)), ...
        qn_canonical_path(intlab_root)), ...
        'startintlab does not resolve inside the requested INTLAB root.');
end
startintlab;
assert(exist('intval','class')==8 || exist('intval','file')==2, ...
    'INTLAB is unavailable.');
assert(exist('veigs','file')==2 && exist('veig','file')==2, ...
    'veigs/veig is unavailable.');
expected_veigs_commit='6556d39a0d9819bb172d232062b698aa76e420f6';
veigs_file=which('veigs'); veig_file=which('veig');
[actual_veigs_commit,veigs_dirty]=qn_git_state(fileparts(veigs_file));
[actual_veig_commit,veig_dirty]=qn_git_state(fileparts(veig_file));
assert(strcmp(actual_veigs_commit,expected_veigs_commit), ...
    'veigs commit mismatch: expected %s, found %s.', ...
    expected_veigs_commit,actual_veigs_commit);
assert(strcmp(actual_veig_commit,expected_veigs_commit), ...
    'veig resolves outside the pinned veigs checkout.');
assert(~veigs_dirty && ~veig_dirty, ...
    'The executed veigs checkout has modifications or untracked files.');
setappdata(0,'qn_veigs_commit',actual_veigs_commit);
setappdata(0,'qn_intlab_start_sha256',qn_sha256_file(intlab_start_file));
setappdata(0,'qn_intlab_tree_sha256', ...
    qn_directory_sha256(fileparts(intlab_start_file)));
end
