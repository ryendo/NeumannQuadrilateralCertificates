function metadata = qn_runtime_provenance()
% Provenance fields embedded in current-schema certificate summaries.

repo_root=fileparts(fileparts(mfilename('fullpath')));
source_commit=strtrim(getenv('QN_SOURCE_COMMIT'));
source_dirty_text=strtrim(getenv('QN_SOURCE_DIRTY'));
[actual_source_commit,source_dirty]=qn_git_state(repo_root,{'results/'});
if isempty(source_commit)
    source_commit=actual_source_commit;
else
    assert(~isempty(regexp(source_commit,'^[0-9a-f]{40}$','once')), ...
        'QN_SOURCE_COMMIT must be a full lowercase Git commit.');
    assert(strcmp(source_commit,actual_source_commit), ...
        'QN_SOURCE_COMMIT does not match the executing checkout.');
    assert(any(strcmp(source_dirty_text,{'0','1'})), ...
        'QN_SOURCE_DIRTY must be 0 or 1 when QN_SOURCE_COMMIT is set.');
    assert(strcmp(source_dirty_text,num2str(double(source_dirty))), ...
        'QN_SOURCE_DIRTY does not match the executing checkout.');
end

expected_veigs_commit='6556d39a0d9819bb172d232062b698aa76e420f6';
if isappdata(0,'qn_veigs_commit')
    veigs_commit=getappdata(0,'qn_veigs_commit');
else
    veigs_file=which('veigs');
    assert(~isempty(veigs_file),'veigs is unavailable for provenance checking.');
    [veigs_commit,veigs_dirty]=qn_git_state(fileparts(veigs_file));
    assert(~veigs_dirty,'The veigs checkout has tracked modifications.');
end
assert(strcmp(veigs_commit,expected_veigs_commit), ...
    'veigs commit mismatch: expected %s, found %s.', ...
    expected_veigs_commit,veigs_commit);

run_id=strtrim(getenv('QN_RUN_ID'));
if isempty(run_id)
    run_id=['interactive-' source_commit(1:12)];
end
assert(~isempty(regexp(run_id,'^[A-Za-z0-9._-]+$','once')), ...
    'QN_RUN_ID contains unsupported characters.');

metadata=struct('source_commit',source_commit, ...
    'source_dirty',logical(source_dirty), ...
    'source_code_sha256',qn_source_code_sha256(repo_root), ...
    'veigs_commit',veigs_commit,'run_id',run_id, ...
    'matlab_version',version,'matlab_release',version('-release'), ...
    'intlab_start_sha256',getappdata(0,'qn_intlab_start_sha256'), ...
    'intlab_tree_sha256',getappdata(0,'qn_intlab_tree_sha256'));
end
