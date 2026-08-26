function repo_root=qn_setup(intlab_root,veigs_root)
% Add the proof code and its two verified numerical dependencies.

repo_root=fileparts(fileparts(mfilename('fullpath')));
if nargin>=1 && ~isempty(intlab_root), addpath(intlab_root); end
if nargin>=2 && ~isempty(veigs_root), addpath(veigs_root); end
addpath(fullfile(repo_root,'src'),fullfile(repo_root,'tests'));
if exist('startintlab','file')==2, startintlab; end
assert(exist('intval','class')==8 || exist('intval','file')==2, ...
    'INTLAB is unavailable.');
assert(exist('veigs','file')==2 && exist('veig','file')==2, ...
    'veigs/veig is unavailable.');
end
