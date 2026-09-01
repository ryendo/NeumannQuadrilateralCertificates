function summary = qn_summarize_local_results(results_dir)
% Summarize the saved output for Algorithm 1 and Table 3.

if nargin<1 || isempty(results_dir)
    root=fileparts(fileparts(mfilename('fullpath')));
    results_dir=fullfile(root,'results','local');
end
files=dir(fullfile(results_dir,'res_*.csv'));
assert(~isempty(files),'No local result files found.');
T=table();
for k=1:numel(files)
    Tk=readtable(fullfile(files(k).folder,files(k).name));
    if k==1, T=Tk; else, T=[T;Tk]; end %#ok<AGROW>
end
assert(height(T)>0,'Local result files contain no data rows.');
ids=sort(T.box_id(:));
unique_ids=unique(ids);
id_coverage_complete=numel(unique_ids)==numel(ids) && ...
    isequal(unique_ids,(1:max(ids))');
worker_ids=zeros(numel(files),1);
for k=1:numel(files)
    token=regexp(files(k).name,'^res_(\d+)\.csv$','tokens','once');
    assert(~isempty(token),'Unexpected local result filename: %s',files(k).name);
    worker_ids(k)=str2double(token{1});
end
done=dir(fullfile(results_dir,'done_*.txt')); done_ids=zeros(numel(done),1);
for k=1:numel(done)
    token=regexp(done(k).name,'^done_(\d+)\.txt$','tokens','once');
    if ~isempty(token), done_ids(k)=str2double(token{1}); end
end
done_complete=isequal(sort(worker_ids),sort(done_ids)) && ...
    numel(unique(worker_ids))==numel(worker_ids);
[minS,iS]=min(T.inf_S); [minL1,iL1]=min(T.lambda1_lower);
[maxL2,iL2]=max(T.lambda2_upper); maxDepth=max(T.max_subdivision_depth);
legacy_verified=all(T.ok==1) && minS>0 && id_coverage_complete && done_complete;
table3_fields={'mass_matrix_lower','lambda3_lower','lambda3_minus_3pi2_over2_lower'};
table3_fields_present=all(ismember(table3_fields,T.Properties.VariableNames));
minMass=nan; minMassBox=nan; minL3=nan; minL3Box=nan; minL3Gap=nan; minL3GapBox=nan;
if table3_fields_present
    [minMass,iMass]=min(T.mass_matrix_lower); minMassBox=T.box_id(iMass);
    [minL3,iL3]=min(T.lambda3_lower); minL3Box=T.box_id(iL3);
    [minL3Gap,iL3Gap]=min(T.lambda3_minus_3pi2_over2_lower);
    minL3GapBox=T.box_id(iL3Gap);
end
verified=legacy_verified && table3_fields_present && minL1>0 && ...
    minMass>0 && minL3>16 && minL3Gap>0;
summary=struct('verified',verified, ...
    'legacy_certificate_verified',legacy_verified, ...
    'table3_fields_present',table3_fields_present, ...
    'result_files',numel(files),'top_level_boxes',height(T), ...
    'certified_boxes',sum(T.ok==1),'failures',sum(T.ok~=1), ...
    'id_coverage_complete',id_coverage_complete,'done_markers_complete',done_complete, ...
    'minimum_certified_lower_bound_for_S',minS,'minimum_S_box_id',T.box_id(iS), ...
    'maximum_subdivision_depth',maxDepth, ...
    'lambda1_lower_min',minL1,'lambda1_lower_min_box_id',T.box_id(iL1), ...
    'lambda2_upper_max',maxL2,'lambda2_upper_max_box_id',T.box_id(iL2), ...
    'mass_matrix_lower_min',minMass,'mass_matrix_lower_min_box_id',minMassBox, ...
    'lambda3_lower_min',minL3,'lambda3_lower_min_box_id',minL3Box, ...
    'lambda3_minus_3pi2_over2_lower_min',minL3Gap, ...
    'lambda3_gap_lower_min_box_id',minL3GapBox, ...
    'wall_time_seconds',max(T.elapsed_seconds));
if summary.verified
    fprintf(['SAVED RESULT CHECK PASSED: all rows pass their recorded checks; ' ...
        'the saved bounds have S(t,e)>0, lambda_1>0, M>0, and lambda_3>16.\n']);
elseif summary.legacy_certificate_verified && ~summary.table3_fields_present
    fprintf(['STALE RESULT SCHEMA: the saved cover verifies the previous S/index-1 ' ...
        'certificate but predates the mass/index-3 output fields.\n']);
else
    fprintf('SAVED RESULT CHECK FAILED: failures=%d, min S=%.17g.\n', ...
        summary.failures,minS);
end
fprintf(['files=%d boxes=%d minS=%.17g depth=%d lambda1>=%.17g ' ...
    'lambda2<=%.17g mass>=%.17g lambda3>=%.17g\n'], ...
    summary.result_files,summary.top_level_boxes,minS,maxDepth,minL1,maxL2,minMass,minL3);
end
