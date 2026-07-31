function summary = qn_summarize_local_results(results_dir)
% Summarize and verify saved DQ2 CSV files using MATLAB only.

if nargin<1 || isempty(results_dir)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
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
[minMass,iMass]=min(T.mass_eigenvalue_lower);
[minL3,iL3]=min(T.lambda3_lower);
[minGap3,iGap3]=min(T.lambda3_gap_lower);
verified=all(T.ok==1) && minS>0 && minMass>0 && minGap3>0 && ...
    id_coverage_complete && done_complete;
summary=struct('verified',verified, ...
    'result_files',numel(files),'top_level_boxes',height(T), ...
    'certified_boxes',sum(T.ok==1),'failures',sum(T.ok~=1), ...
    'id_coverage_complete',id_coverage_complete,'done_markers_complete',done_complete, ...
    'minimum_certified_lower_bound_for_S',minS,'minimum_S_box_id',T.box_id(iS), ...
    'maximum_subdivision_depth',maxDepth, ...
    'lambda1_lower_min',minL1,'lambda1_lower_min_box_id',T.box_id(iL1), ...
    'lambda2_upper_max',maxL2,'lambda2_upper_max_box_id',T.box_id(iL2), ...
    'mass_eigenvalue_lower_min',minMass,'mass_eigenvalue_lower_min_box_id',T.box_id(iMass), ...
    'lambda3_lower_min',minL3,'lambda3_lower_min_box_id',T.box_id(iL3), ...
    'lambda3_gap_lower_min',minGap3,'lambda3_gap_lower_min_box_id',T.box_id(iGap3), ...
    'wall_time_seconds',max(T.elapsed_seconds));
if summary.verified
    fprintf('VERIFIED: S(t,e)>0 is certified on all saved boxes.\n');
else
    fprintf('NOT VERIFIED: failures=%d, min S=%.17g.\n',summary.failures,minS);
end
fprintf(['files=%d boxes=%d minS=%.17g depth=%d lambda1>=%.17g ' ...
    'lambda2<=%.17g lambda_min(M)>=%.17g lambda3>=%.17g gap3>=%.17g\n'], ...
    summary.result_files,summary.top_level_boxes,minS,maxDepth,minL1,maxL2, ...
    minMass,minL3,minGap3);
end
