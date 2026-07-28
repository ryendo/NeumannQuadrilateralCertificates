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
% Preserve a genuine column stack after linear worker slicing.  If a row cell
% reaches a two-subscript append, MATLAB can expand the other columns with
% empty cells; a later pop would then return [] instead of a box struct.
stack=stack(:);
initial_retained=numel(stack);
cert=0; disc=0; bis=0; unv=0; two=0; veigs_cert=0; max_d=0;
min_margin=inf; min_at=[]; prior_wall=0;
unverified_boxes=struct('center',{},'half_widths',{},'depth',{}, ...
    'reason',{},'gap',{},'split_dim',{});
checkpoint_file='';
if ~isempty(output_file)
    checkpoint_file=[output_file '.checkpoint.mat'];
    if isfile(checkpoint_file)
        saved=load(checkpoint_file,'state'); state=saved.state;
        assert(state.version==2 && state.worker_id==worker_id && ...
            state.worker_count==worker_count && state.n_init==n_init && ...
            state.max_depth==max_depth,'qn:CheckpointMismatch', ...
            'Global-cover checkpoint does not match this worker invocation.');
        stack=state.stack(:); cert=state.cert; disc=state.disc;
        bis=state.bis; unv=state.unv; two=state.two; max_d=state.max_d;
        if isfield(state,'veigs_cert'), veigs_cert=state.veigs_cert; end
        min_margin=state.min_margin; min_at=state.min_at;
        unverified_boxes=state.unverified_boxes; prior_wall=state.wall_seconds;
        fprintf('RESUME global worker=%d cert=%d disc=%d bis=%d stack=%d wall_s=%.1f\n', ...
            worker_id,cert,disc,bis,numel(stack),prior_wall);
    end
end
t0=tic;
while ~isempty(stack)
    B=stack{end,1}; stack(end,:)=[];
    assert(isstruct(B) && isscalar(B) && isfield(B,'depth'), ...
        'qn:StackCorrupt','Global cover stack contains a non-box entry.');
    max_d=max(max_d,B.depth);
    if qn_box_inside_ball(B.center,B.half_widths,C.rho_seam) || ...
            qn_box_outside_pk(B.center,B.half_widths)
        disc=disc+1;
        maybe_checkpoint(checkpoint_file,stack,cert,disc,bis,unv,two,max_d, ...
            veigs_cert,min_margin,min_at,unverified_boxes,prior_wall+toc(t0), ...
            worker_id,worker_count,n_init,max_depth);
        continue;
    end
    [verdict,info,gap,sdim]=qn_certify_box(B.center,B.half_widths);
    if strcmp(verdict,'cert')
        cert=cert+1; veigs_cert=veigs_cert+strcmp(info.route,'veigs');
        if info.margin<min_margin, min_margin=info.margin; min_at=B.center'; end
    elseif B.depth<max_depth
        [L,R]=qn_bisect_box(B,sdim); stack{end+1,1}=L; stack{end+1,1}=R; bis=bis+1;
    else
        unv=unv+1;
        if isfield(info,'reason'), reason=info.reason; else, reason='unknown'; end
        unverified_boxes(end+1)=struct('center',B.center', ... %#ok<AGROW>
            'half_widths',B.half_widths','depth',B.depth,'reason',reason, ...
            'gap',gap,'split_dim',sdim);
        fprintf(['UNVERIFIED depth=%d reason=%s center=[%.17g %.17g %.17g %.17g] ' ...
            'half_widths=[%.17g %.17g %.17g %.17g]\n'],B.depth,reason, ...
            B.center(1),B.center(2),B.center(3),B.center(4), ...
            B.half_widths(1),B.half_widths(2),B.half_widths(3),B.half_widths(4));
    end
    if verbose && mod(cert+disc+bis,200)==0
        fprintf('global: cert=%d disc=%d bis=%d unv=%d stack=%d\n',cert,disc,bis,unv,numel(stack));
    end
    maybe_checkpoint(checkpoint_file,stack,cert,disc,bis,unv,two,max_d, ...
        veigs_cert,min_margin,min_at,unverified_boxes,prior_wall+toc(t0), ...
        worker_id,worker_count,n_init,max_depth);
end
result=struct('verified',cert,'discarded',disc,'bisected',bis, ...
    'unverified',unv,'two_vector',two,'veigs_certified',veigs_cert,'max_depth',max_d, ...
    'min_certified_margin',min_margin,'min_margin_at',min_at, ...
    'wall_seconds',prior_wall+toc(t0),'n_init',n_init,'initial_retained',initial_retained, ...
    'initial_retained_all',initial_retained_all,'worker_id',worker_id, ...
    'worker_count',worker_count,'rho_seam',mid(C.rho_seam), ...
    'complete',unv==0 && veigs_cert==cert);
result.unverified_boxes=unverified_boxes;
fprintf(['RESULT global verified=%d discarded=%d bisected=%d unverified=%d ' ...
    'two_vector=%d max_depth=%d min_certified_margin=%.17g wall_s=%.1f\n'], ...
    cert,disc,bis,unv,two,max_d,min_margin,result.wall_seconds);
if ~isempty(output_file)
    temporary_output=[output_file '.tmp'];
    fid=fopen(temporary_output,'w'); assert(fid>=0,'Cannot open output file.');
    cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));
    clear cleanup
    movefile(temporary_output,output_file,'f');
    if isfile(checkpoint_file), delete(checkpoint_file); end
end
end

function maybe_checkpoint(filename,stack,cert,disc,bis,unv,two,max_d, ...
        veigs_cert,min_margin,min_at,unverified_boxes,wall_seconds,worker_id,worker_count,n_init,max_depth)
if isempty(filename) || mod(cert+disc+bis+unv,10)~=0, return; end
state=struct('version',2,'stack',{stack(:)},'cert',cert,'disc',disc, ...
    'bis',bis,'unv',unv,'two',two,'veigs_cert',veigs_cert, ...
    'max_d',max_d,'min_margin',min_margin, ...
    'min_at',min_at,'unverified_boxes',unverified_boxes, ...
    'wall_seconds',wall_seconds,'worker_id',worker_id, ...
    'worker_count',worker_count,'n_init',n_init,'max_depth',max_depth);
temporary=[filename '.tmp.mat'];
save(temporary,'state','-v7');
movefile(temporary,filename,'f');
end
