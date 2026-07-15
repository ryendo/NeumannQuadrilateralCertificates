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
[minS,iS]=min(T.inf_S); [minL1,iL1]=min(T.lambda1_lower);
[maxL2,iL2]=max(T.lambda2_upper); maxDepth=max(T.max_subdivision_depth);
summary=struct('verified',all(T.ok==1)&&minS>0, ...
    'result_files',numel(files),'top_level_boxes',height(T), ...
    'certified_boxes',sum(T.ok==1),'failures',sum(T.ok~=1), ...
    'minimum_certified_lower_bound_for_S',minS,'minimum_S_box_id',T.box_id(iS), ...
    'maximum_subdivision_depth',maxDepth, ...
    'lambda1_lower_min',minL1,'lambda1_lower_min_box_id',T.box_id(iL1), ...
    'lambda2_upper_max',maxL2,'lambda2_upper_max_box_id',T.box_id(iL2), ...
    'wall_time_seconds',max(T.elapsed_seconds));
if summary.verified
    fprintf('VERIFIED: S(t,e)>0 is certified on all saved boxes.\n');
else
    fprintf('NOT VERIFIED: failures=%d, min S=%.17g.\n',summary.failures,minS);
end
fprintf('files=%d boxes=%d minS=%.17g depth=%d lambda1>=%.17g lambda2<=%.17g\n', ...
    summary.result_files,summary.top_level_boxes,minS,maxDepth,minL1,maxL2);
end
