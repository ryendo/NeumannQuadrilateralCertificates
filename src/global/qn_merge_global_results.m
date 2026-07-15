function result = qn_merge_global_results(results_dir,worker_count,output_file)
% Merge independent root-box worker summaries. No proof arithmetic is redone.

if nargin<2, worker_count=10; end
if nargin<3, output_file=fullfile(results_dir,'summary.json'); end
items=cell(worker_count,1);
for k=1:worker_count
    path=fullfile(results_dir,sprintf('worker_%03d.json',k));
    assert(isfile(path),'Missing worker result: %s',path);
    items{k}=jsondecode(fileread(path));
    assert(items{k}.worker_id==k,'Unexpected worker_id in %s.',path);
    assert(items{k}.worker_count==worker_count,'Worker-count mismatch in %s.',path);
end
assert(all(cellfun(@(x)x.n_init==items{1}.n_init,items)),'n_init differs across workers.');
assert(all(cellfun(@(x)x.initial_retained_all==items{1}.initial_retained_all,items)), ...
    'Initial-cover size differs across workers.');
fields={'verified','discarded','bisected','unverified','two_vector','initial_retained'};
result=struct();
for f=1:numel(fields)
    name=fields{f}; total=0;
    for k=1:worker_count, total=total+items{k}.(name); end
    result.(name)=total;
end
result.max_depth=max(cellfun(@(x)x.max_depth,items));
margins=inf(worker_count,1);
for k=1:worker_count
    if items{k}.verified>0 && isnumeric(items{k}.min_certified_margin) && ...
            isscalar(items{k}.min_certified_margin) && isfinite(items{k}.min_certified_margin)
        margins(k)=items{k}.min_certified_margin;
    end
end
[result.min_certified_margin,kmin]=min(margins);
if isfinite(result.min_certified_margin)
    result.min_margin_at=items{kmin}.min_margin_at;
else
    result.min_margin_at=[];
end
result.wall_seconds_max=max(cellfun(@(x)x.wall_seconds,items));
result.worker_count=worker_count;
result.n_init=items{1}.n_init;
result.initial_retained_all=items{1}.initial_retained_all;
result.rho_seam=items{1}.rho_seam;
all_workers_complete=all(cellfun(@(x)logical(x.complete),items));
result.complete=all_workers_complete && result.unverified==0 && ...
    result.initial_retained==result.initial_retained_all && isfinite(result.min_certified_margin);
fid=fopen(output_file,'w'); assert(fid>=0,'Cannot open merged output file.');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));
fprintf(['MERGED global verified=%d discarded=%d bisected=%d unverified=%d ' ...
    'max_depth=%d min_margin=%.17g complete=%d\n'],result.verified,result.discarded, ...
    result.bisected,result.unverified,result.max_depth,result.min_certified_margin,result.complete);
end
