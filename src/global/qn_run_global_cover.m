function result = qn_run_global_cover(n_init,max_depth,verbose,output_file,worker_id,worker_count)
% Serial certified adaptive cover for the paper's global step.

if nargin<1 || isempty(n_init), n_init=3; end
C=qn_global_constants();
if nargin<2 || isempty(max_depth), max_depth=C.max_depth; end
if nargin<3 || isempty(verbose), verbose=true; end
if nargin<4, output_file=''; end
if nargin<5 || isempty(worker_id), worker_id=1; end
if nargin<6 || isempty(worker_count), worker_count=1; end

allboxes=qn_initial_cover(n_init); stack=cell(0,1);
for k=1:numel(allboxes)
    B=allboxes{k};
    if qn_d4_is_canonical(B.center) && ~qn_box_outside_pk(B.center,B.half_widths)
        stack{end+1,1}=B; %#ok<AGROW>
    end
end
initial_retained_all=numel(stack);
stack=stack(worker_id:worker_count:end);
initial_retained=numel(stack);
cert=0; disc=0; bis=0; unv=0; two=0; max_d=0;
min_margin=inf; min_at=[]; t0=tic;
while ~isempty(stack)
    B=stack{end}; stack(end)=[]; max_d=max(max_d,B.depth);
    if qn_box_inside_ball(B.center,B.half_widths,C.rho_seam) || ...
            qn_box_outside_pk(B.center,B.half_widths)
        disc=disc+1; continue;
    end
    [verdict,info,~,sdim]=qn_certify_box(B.center,B.half_widths);
    if strcmp(verdict,'cert')
        cert=cert+1; two=two+strcmp(info.route,'2v');
        if info.margin<min_margin, min_margin=info.margin; min_at=B.center'; end
    elseif B.depth<max_depth
        [L,R]=qn_bisect_box(B,sdim); stack{end+1,1}=L; stack{end+1,1}=R; bis=bis+1;
    else
        unv=unv+1;
    end
    if verbose && mod(cert+disc+bis,200)==0
        fprintf('global: cert=%d disc=%d bis=%d unv=%d stack=%d\n',cert,disc,bis,unv,numel(stack));
    end
end
result=struct('verified',cert,'discarded',disc,'bisected',bis, ...
    'unverified',unv,'two_vector',two,'max_depth',max_d, ...
    'min_certified_margin',min_margin,'min_margin_at',min_at, ...
    'wall_seconds',toc(t0),'n_init',n_init,'initial_retained',initial_retained, ...
    'initial_retained_all',initial_retained_all,'worker_id',worker_id, ...
    'worker_count',worker_count,'rho_seam',mid(C.rho_seam),'complete',unv==0);
fprintf(['RESULT global verified=%d discarded=%d bisected=%d unverified=%d ' ...
    'two_vector=%d max_depth=%d min_certified_margin=%.17g wall_s=%.1f\n'], ...
    cert,disc,bis,unv,two,max_d,min_margin,result.wall_seconds);
if ~isempty(output_file)
    fid=fopen(output_file,'w'); assert(fid>=0,'Cannot open output file.');
    cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));
end
end
