function result = qn_global_certified_cover(n_init,max_depth,verbose,output_file,worker_id,worker_count)
% Implements the Appendix B finite-cover procedure for Theorem 6.2,
% including the acceptance, discard, subdivision, and completion conditions
% in (44)-(49).

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
initial_box_ids=zeros(initial_retained,1);
for k=1:initial_retained
    initial_box_ids(k)=stack{k}.initial_id;
end
cert=0; disc=0; bis=0; unv=0; max_d=0;
delta_star_lower=inf; delta_star_at=[];
q_lower_min=inf; lambda_lower_min=inf;
unverified_boxes=struct('center',{},'half_widths',{},'depth',{}, ...
    'reason',{},'split_dim',{});
t0=tic;
while ~isempty(stack)
    B=stack{end,1}; stack(end,:)=[];
    assert(isstruct(B) && isscalar(B) && isfield(B,'depth'), ...
        'qn:StackCorrupt','Global cover stack contains a non-box entry.');
    max_d=max(max_d,B.depth);
    if qn_box_inside_ball(B.center,B.half_widths,C.rho_sharp_over_2) || ...
            qn_box_outside_pk(B.center,B.half_widths)
        disc=disc+1;
        continue;
    end
    [verdict,info,sdim]=qn_certify_box(B.center,B.half_widths);
    if strcmp(verdict,'cert')
        cert=cert+1;
        q_lower_min=min(q_lower_min,info.q_lower);
        lambda_lower_min=min(lambda_lower_min,info.lambda1(1));
        if info.delta_lower<delta_star_lower
            delta_star_lower=info.delta_lower;
            delta_star_at=B.center';
        end
    elseif B.depth<max_depth
        [L,R]=qn_bisect_box(B,sdim); stack{end+1,1}=L; stack{end+1,1}=R; bis=bis+1;
    else
        unv=unv+1;
        if isfield(info,'reason'), reason=info.reason; else, reason='unknown'; end
        unverified_boxes(end+1)=struct('center',B.center', ... %#ok<AGROW>
            'half_widths',B.half_widths','depth',B.depth,'reason',reason, ...
            'split_dim',sdim);
        fprintf(['UNVERIFIED depth=%d reason=%s center=[%.17g %.17g %.17g %.17g] ' ...
            'half_widths=[%.17g %.17g %.17g %.17g]\n'],B.depth,reason, ...
            B.center(1),B.center(2),B.center(3),B.center(4), ...
            B.half_widths(1),B.half_widths(2),B.half_widths(3),B.half_widths(4));
    end
    if verbose && mod(cert+disc+bis,200)==0
        fprintf('global: cert=%d disc=%d bis=%d unv=%d stack=%d\n',cert,disc,bis,unv,numel(stack));
    end
end
worker_complete=unv==0;
covers_all_initial=initial_retained==initial_retained_all;
complete=worker_complete && covers_all_initial && ...
    isfinite(q_lower_min) && q_lower_min>0 && ...
    isfinite(lambda_lower_min) && lambda_lower_min>0 && ...
    isfinite(delta_star_lower) && delta_star_lower>0;
result=struct('verified',cert,'discarded',disc,'bisected',bis, ...
    'unverified',unv,'veigs_certified',cert,'max_depth',max_d, ...
    'delta_star_lower',delta_star_lower,'delta_star_at',delta_star_at, ...
    'q_lower_min',q_lower_min,'lambda_lower_min',lambda_lower_min, ...
    'min_certified_margin',delta_star_lower,'min_margin_at',delta_star_at, ...
    'wall_seconds',toc(t0),'n_init',n_init,'initial_retained',initial_retained, ...
    'initial_retained_all',initial_retained_all,'worker_id',worker_id, ...
    'worker_count',worker_count,'rho_sharp_over_2',mid(C.rho_sharp_over_2), ...
    'schema_version',2,'initial_box_ids',initial_box_ids', ...
    'worker_complete',worker_complete,'complete',complete);
result.unverified_boxes=unverified_boxes;
fprintf(['RESULT global verified=%d discarded=%d bisected=%d unverified=%d ' ...
    'max_depth=%d delta_star_lower=%.17g wall_s=%.1f\n'], ...
    cert,disc,bis,unv,max_d,delta_star_lower,result.wall_seconds);
if ~isempty(output_file)
    temporary_output=[output_file '.tmp'];
    fid=fopen(temporary_output,'w'); assert(fid>=0,'Cannot open output file.');
    cleanup=onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));
    clear cleanup
    movefile(temporary_output,output_file,'f');
end
end
