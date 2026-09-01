function result = qn_merge_global_results(results_dir,worker_count,output_file)
% Merge worker summaries and verify the completion conditions in (48).

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
assert(all(cellfun(@(x)x.veigs_certified==x.verified,items)), ...
    'A worker summary has inconsistent accepted and veigs-certified counts.');
reference_radius=reported_radius(items{1});
assert(all(cellfun(@(x)reported_radius(x)==reference_radius,items)), ...
    'The radius rho^sharp/2 differs across workers.');
assert(all(cellfun(@tree_accounting_holds,items)), ...
    'A worker summary violates the binary-tree accounting for (47).');
fields={'verified','discarded','bisected','unverified','veigs_certified', ...
    'initial_retained'};
result=struct();
for f=1:numel(fields)
    name=fields{f};
    result.(name)=sum(cellfun(@(item)item.(name),items));
end
result.max_depth=max(cellfun(@(x)x.max_depth,items));
current_schema_present=all(cellfun(@has_current_schema,items));
all_delta_records_valid=all(cellfun(@has_valid_delta_record,items));
current_worker_bounds_valid=current_schema_present && ...
    all(cellfun(@has_valid_current_bounds,items));
[current_initial_assignment_valid,current_initial_box_ids,expected_initial_box_ids]= ...
    validate_initial_assignments(items,current_schema_present, ...
    items{1}.initial_retained_all);
deltas=inf(worker_count,1);
q_lowers=inf(worker_count,1);
lambda_lowers=inf(worker_count,1);
for k=1:worker_count
    if isfield(items{k},'delta_star_lower')
        value=items{k}.delta_star_lower;
    else
        value=items{k}.min_certified_margin;
    end
    if items{k}.verified>0 && isnumeric(value) && isscalar(value) && isfinite(value)
        deltas(k)=value;
    end
    if current_schema_present && items{k}.verified>0 && ...
            isnumeric(items{k}.q_lower_min) && isscalar(items{k}.q_lower_min) && ...
            isfinite(items{k}.q_lower_min) && ...
            isnumeric(items{k}.lambda_lower_min) && ...
            isscalar(items{k}.lambda_lower_min) && ...
            isfinite(items{k}.lambda_lower_min)
        q_lowers(k)=items{k}.q_lower_min;
        lambda_lowers(k)=items{k}.lambda_lower_min;
    end
end
[result.delta_star_lower,kmin]=min(deltas);
result.min_certified_margin=result.delta_star_lower;
if current_schema_present
    result.q_lower_min=min(q_lowers);
    result.lambda_lower_min=min(lambda_lowers);
else
    result.q_lower_min=[];
    result.lambda_lower_min=[];
end
if isfinite(result.delta_star_lower)
    if isfield(items{kmin},'delta_star_at')
        result.delta_star_at=items{kmin}.delta_star_at;
    else
        result.delta_star_at=items{kmin}.min_margin_at;
    end
    result.min_margin_at=result.delta_star_at;
else
    result.delta_star_at=[];
    result.min_margin_at=[];
end
result.wall_seconds_max=max(cellfun(@(x)x.wall_seconds,items));
result.worker_count=worker_count;
result.n_init=items{1}.n_init;
result.initial_retained_all=items{1}.initial_retained_all;
if isfield(items{1},'rho_sharp_over_2')
    result.rho_sharp_over_2=items{1}.rho_sharp_over_2;
else
    result.rho_sharp_over_2=items{1}.rho_seam;
end
all_workers_complete=all(cellfun(@worker_finished,items));
legacy_completion_recorded=all_workers_complete && result.unverified==0 && ...
    result.initial_retained==result.initial_retained_all && ...
    all_delta_records_valid && ...
    isfinite(result.delta_star_lower) && result.delta_star_lower>0;
result.current_schema_present=current_schema_present;
result.current_worker_bounds_valid=current_worker_bounds_valid;
result.current_initial_assignment_valid=current_initial_assignment_valid;
result.initial_box_ids=current_initial_box_ids';
result.expected_initial_box_ids=expected_initial_box_ids';
if current_schema_present, result.schema_version=2; else, result.schema_version=[]; end
result.legacy_completion_recorded=legacy_completion_recorded;
result.complete=current_schema_present && current_worker_bounds_valid && ...
    current_initial_assignment_valid && ...
    legacy_completion_recorded && ...
    isfinite(result.q_lower_min) && result.q_lower_min>0 && ...
    isfinite(result.lambda_lower_min) && result.lambda_lower_min>0;
fid=fopen(output_file,'w'); assert(fid>=0,'Cannot open merged output file.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));
fprintf(['MERGED global verified=%d discarded=%d bisected=%d unverified=%d ' ...
    'max_depth=%d delta_star_lower=%.17g complete=%d\n'], ...
    result.verified,result.discarded,result.bisected,result.unverified, ...
    result.max_depth,result.delta_star_lower,result.complete);
end

function tf=worker_finished(item)
% New summaries distinguish a finished worker slice from global condition (48).
if isfield(item,'worker_complete')
    tf=logical(item.worker_complete);
else
    tf=logical(item.complete);
end
end

function tf=has_current_schema(item)
required={'delta_star_lower','q_lower_min','lambda_lower_min','worker_complete', ...
    'schema_version','initial_box_ids'};
tf=all(cellfun(@(name)isfield(item,name),required)) && ...
    isnumeric(item.schema_version) && isscalar(item.schema_version) && ...
    item.schema_version==2;
end

function tf=has_valid_delta_record(item)
if item.verified==0
    tf=true;
    return
end
if isfield(item,'delta_star_lower')
    value=item.delta_star_lower;
else
    value=item.min_certified_margin;
end
tf=isnumeric(value) && isscalar(value) && isfinite(value) && value>0;
end

function tf=has_valid_current_bounds(item)
if item.verified==0
    tf=true;
    return
end
tf=isnumeric(item.q_lower_min) && isscalar(item.q_lower_min) && ...
    isfinite(item.q_lower_min) && item.q_lower_min>0 && ...
    isnumeric(item.lambda_lower_min) && isscalar(item.lambda_lower_min) && ...
    isfinite(item.lambda_lower_min) && item.lambda_lower_min>0 && ...
    has_valid_delta_record(item);
end

function tf=tree_accounting_holds(item)
% Each bisection replaces one active box by two children.
tf=item.verified+item.discarded+item.unverified == ...
    item.initial_retained+item.bisected;
end

function value=reported_radius(item)
if isfield(item,'rho_sharp_over_2')
    value=item.rho_sharp_over_2;
else
    value=item.rho_seam; % Field name in the recorded earlier JSON schema.
end
end

function [tf,assigned,expected]=validate_initial_assignments( ...
    items,current_fields_present,total)
assigned=[]; expected=[];
if ~current_fields_present
    tf=false;
    return
end
allboxes=qn_initial_cover(items{1}.n_init);
for k=1:numel(allboxes)
    box=allboxes{k};
    if qn_d4_is_canonical(box.center) && ...
            ~qn_box_outside_pk(box.center,box.half_widths)
        expected(end+1,1)=box.initial_id; %#ok<AGROW>
    end
end
tf=true;
for k=1:numel(items)
    ids=items{k}.initial_box_ids;
    if ~isnumeric(ids) || any(~isfinite(ids(:))) || ...
            any(ids(:)~=floor(ids(:))) || any(ids(:)<1) || ...
            any(ids(:)>items{k}.n_init^4) || ...
            numel(ids)~=items{k}.initial_retained
        tf=false;
        return
    end
    assigned=[assigned;ids(:)]; %#ok<AGROW>
end
tf=total==numel(expected) && numel(assigned)==total && ...
    numel(unique(assigned))==total && ...
    isequal(sort(assigned),sort(expected));
end
