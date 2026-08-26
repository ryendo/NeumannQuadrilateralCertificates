function result = qn_global_root_subcover(root_id,presplit_depth,subworker_id,subworker_count,max_depth,output_file)
% Proof-equivalent parallel subdivision of one retained n_init=3 root box.
% This temporary HPC driver reuses the repository's certified box predicates
% and longest-side bisection.  The 2^presplit_depth leaves form an exact cover.

if nargin<5 || isempty(max_depth), max_depth=60; end
if nargin<6, output_file=''; end
C=qn_global_constants();
allboxes=qn_initial_cover(3); roots=cell(0,1);
for k=1:numel(allboxes)
    B=allboxes{k};
    if qn_d4_is_canonical(B.center) && ~qn_box_outside_pk(B.center,B.half_widths)
        roots{end+1,1}=B; %#ok<AGROW>
    end
end
assert(numel(roots)==16 && root_id>=1 && root_id<=16,'Unexpected retained-root family.');

leaves=roots(root_id);
for level=1:presplit_depth
    next=cell(2*numel(leaves),1);
    for k=1:numel(leaves)
        B=leaves{k}; [~,split_dim]=max(B.half_widths);
        [L,R]=qn_bisect_box(B,split_dim);
        next{2*k-1}=L; next{2*k}=R;
    end
    leaves=next;
end
assert(subworker_id>=1 && subworker_id<=subworker_count,'Invalid subworker ID.');
stack=leaves(subworker_id:subworker_count:end); stack=stack(:);
assigned_leaves=numel(stack);

cert=0; disc=0; bis=0; unv=0; max_d=0; min_margin=inf; min_at=[];
unverified_boxes=struct('center',{},'half_widths',{},'depth',{},'reason',{},'split_dim',{});
t0=tic;
while ~isempty(stack)
    B=stack{end,1}; stack(end,:)=[]; max_d=max(max_d,B.depth);
    if qn_box_inside_ball(B.center,B.half_widths,C.rho_seam) || ...
            qn_box_outside_pk(B.center,B.half_widths)
        disc=disc+1; continue;
    end
    [verdict,info,split_dim]=qn_certify_box(B.center,B.half_widths);
    if strcmp(verdict,'cert')
        cert=cert+1;
        if info.margin<min_margin, min_margin=info.margin; min_at=B.center'; end
    elseif B.depth<max_depth
        [L,R]=qn_bisect_box(B,split_dim);
        stack{end+1,1}=L; stack{end+1,1}=R; bis=bis+1;
    else
        unv=unv+1;
        if isfield(info,'reason'), reason=info.reason; else, reason='unknown'; end
        unverified_boxes(end+1)=struct('center',B.center','half_widths',B.half_widths', ... %#ok<AGROW>
            'depth',B.depth,'reason',reason,'split_dim',split_dim);
    end
    if mod(cert+disc+bis,200)==0
        fprintf('root=%d sub=%d cert=%d disc=%d bis=%d unv=%d stack=%d\n', ...
            root_id,subworker_id,cert,disc,bis,unv,numel(stack));
    end
end
result=struct('verified',cert,'discarded',disc,'bisected',bis,'unverified',unv, ...
    'veigs_certified',cert,'max_depth',max_d,'min_certified_margin',min_margin, ...
    'min_margin_at',min_at,'wall_seconds',toc(t0),'n_init',3, ...
    'initial_retained',1,'initial_retained_all',16,'worker_id',root_id, ...
    'worker_count',18,'rho_seam',mid(C.rho_seam),'complete',unv==0, ...
    'parallel_presplit_depth',presplit_depth,'parallel_subworker_id',subworker_id, ...
    'parallel_subworker_count',subworker_count,'assigned_presplit_leaves',assigned_leaves, ...
    'unverified_boxes',unverified_boxes);
fprintf('RESULT root=%d sub=%d verified=%d discarded=%d bisected=%d unverified=%d max_depth=%d margin=%.17g\n', ...
    root_id,subworker_id,cert,disc,bis,unv,max_d,min_margin);
if ~isempty(output_file)
    temporary_output=[output_file '.tmp'];
    fid=fopen(temporary_output,'w'); assert(fid>=0,'Cannot open output file.');
    cleanup=onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));
    clear cleanup
    movefile(temporary_output,output_file,'f');
end
end
