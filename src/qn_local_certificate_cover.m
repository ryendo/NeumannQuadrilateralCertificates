function qn_local_certificate_cover(worker_id, worker_count, face_subdivisions, outdir)
% Process one worker slice of the local Algorithm 1 certificate cover.
% This is the covering/subdivision loop described in Appendix A.4.
%
% PDF notation:
%   local target f_2(p)>pi^(-2); p=t e                  (21), Appendix A.1
%   direction charts and radial intervals               (40)-(41)
%   S(t,e)>0: definition and target                     (23), (26)
%   Single-box certificate returns VERIFIED/UNDECIDED
%
% Evaluation:
%   The face boxes E_j are ordinary floating-point endpoints, immediately
%   converted to INTLAB intervals inside qn_single_box_certificate.m. Every
%   undecided product box is subdivided and the complete certificate is rerun.
%
% Writes per-box result lines to outdir/res_<procid>.csv and a DONE marker.
% Every real-valued CSV field is written with %.17e so the decimal strings
% retain round-trip double precision for the INTLAB interval endpoints.

if nargin < 4, outdir = 'results'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
provenance=qn_runtime_provenance();
PI = intval('pi');
rho = 3232/(27*PI^6);                         % rho# in Sec. 1
tmax = intval(sup(rho))*intval('1.000000000001');
radial_grid = tmax * (intval([0 2 4 6 8 10 12 14 15 16])/intval(16)); % T intervals, refined near ρ#
result_path=sprintf('%s/res_%03d.csv',outdir,worker_id);
done_path=sprintf('%s/done_%03d.txt',outdir,worker_id);
fid = fopen(result_path, 'w');
fprintf(fid, ['box_id,ok,inf_S,max_subdivision_depth,elapsed_seconds,' ...
    'lambda1_lower,lambda2_upper,mass_matrix_lower,lambda3_lower,' ...
    'lambda3_minus_3pi2_over2_lower\n']);
face_box_table = face_boxes(face_subdivisions); % E boxes for B in Eq. (40)
total_boxes = size(face_box_table, 1);
tstart = tic;
failed_boxes = 0; worst_lower_S = inf; worst_mass_lower = inf;
worst_lambda3_lower = inf; worst_lambda3_gap_lower = inf;
for box_id = worker_id:worker_count:total_boxes
    face_box = face_box_table(box_id, :);       % one E; radial_grid supplies T
    proof = certify_with_subdivision(face_box,0,radial_grid,1:(size(radial_grid,2)-1));
    worst_lower_S = min(worst_lower_S, proof.lower_S);
    worst_mass_lower = min(worst_mass_lower,proof.mass_lower);
    worst_lambda3_lower = min(worst_lambda3_lower,proof.lambda3_lower);
    worst_lambda3_gap_lower = min(worst_lambda3_gap_lower,proof.lambda3_gap_lower);
    if ~proof.verified, failed_boxes = failed_boxes + 1; end
    fprintf(fid, '%d,%d,%.17e,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n', ...
        box_id, proof.verified, proof.lower_S, proof.depth, toc(tstart), ...
        proof.lambda1_lower, proof.lambda2_upper, proof.mass_lower, ...
        proof.lambda3_lower, proof.lambda3_gap_lower);
    if mod(box_id, 40*worker_count) < worker_count
        fprintf('proc %d: box %d/%d worstS=%.3f fails=%d (%.0fs)\n', ...
                worker_id, box_id, total_boxes, worst_lower_S, failed_boxes, toc(tstart));
    end
end
fclose(fid);
fid2 = fopen(done_path, 'w');
fprintf(fid2, ['worstS=%.17e worstMass=%.17e worstLambda3=%.17e ' ...
    'worstLambda3Gap=%.17e nfail=%d time=%.17e\n'],worst_lower_S, ...
    worst_mass_lower,worst_lambda3_lower,worst_lambda3_gap_lower, ...
    failed_boxes,toc(tstart));
fclose(fid2);
metadata=provenance;
metadata.schema_version=1;
metadata.worker_id=worker_id;
metadata.worker_count=worker_count;
metadata.face_subdivisions=face_subdivisions;
metadata.expected_top_level_boxes=total_boxes;
metadata.result_sha256=qn_sha256_file(result_path);
metadata.done_sha256=qn_sha256_file(done_path);
metadata_path=sprintf('%s/meta_%03d.json',outdir,worker_id);
temporary_metadata=[metadata_path '.tmp'];
fid3=fopen(temporary_metadata,'w'); assert(fid3>=0,'Cannot write local metadata.');
cleanup=onCleanup(@() fclose(fid3));
fprintf(fid3,'%s\n',jsonencode(metadata,'PrettyPrint',true));
clear cleanup
movefile(temporary_metadata,metadata_path,'f');
end

function proof = certify_with_subdivision(face_box,depth,radial_grid,t_levels)
maxdepth = 6;
res = qn_single_box_certificate(face_box,radial_grid,t_levels);
proof = summarize_results(res,depth,t_levels);
if proof.verified || depth >= maxdepth, return; end
failed = t_levels(~[res(t_levels).ok]);
proof.verified = true;
for child_box = subdivide_face_box(face_box)'
    child = certify_with_subdivision(child_box',depth+1,radial_grid,failed);
    proof = merge_proofs(proof,child);
    if ~proof.verified && depth >= maxdepth-1, return; end
end
end

function proof = summarize_results(res,depth,t_levels)
ok=[res(t_levels).ok];
accepted=t_levels(ok);
lower=min([res(accepted).infS]);
if isempty(lower), lower=inf; end
proof=struct('verified',all(ok),'lower_S',lower,'depth',depth, ...
    'lambda1_lower',inf,'lambda2_upper',-inf,'mass_lower',inf, ...
    'lambda3_lower',inf,'lambda3_gap_lower',inf);
if ~isfield(res,'lam1'), return; end
for k=accepted
    if ~isempty(res(k).lam1)
        proof.lambda1_lower=min(proof.lambda1_lower,res(k).lam1(1));
        proof.lambda2_upper=max(proof.lambda2_upper,res(k).lam2(2));
    end
    if res(k).ok && isfield(res,'lam3') && ~isempty(res(k).lam3)
        proof.mass_lower=min(proof.mass_lower,res(k).mass_lower);
        proof.lambda3_lower=min(proof.lambda3_lower,res(k).lam3(1));
        proof.lambda3_gap_lower=min(proof.lambda3_gap_lower,res(k).lambda3_gap_lower);
    end
end
end

function proof = merge_proofs(proof,child)
proof.verified=proof.verified && child.verified;
proof.lower_S=min(proof.lower_S,child.lower_S);
proof.depth=max(proof.depth,child.depth);
proof.lambda1_lower=min(proof.lambda1_lower,child.lambda1_lower);
proof.lambda2_upper=max(proof.lambda2_upper,child.lambda2_upper);
proof.mass_lower=min(proof.mass_lower,child.mass_lower);
proof.lambda3_lower=min(proof.lambda3_lower,child.lambda3_lower);
proof.lambda3_gap_lower=min(proof.lambda3_gap_lower,child.lambda3_gap_lower);
end

function child_boxes = subdivide_face_box(parent_box)
lo = parent_box(3:2:7); hi = parent_box(4:2:8);
md = (lo + hi)/2;
child_boxes = zeros(8, 8);
for c = 0:7
    ch = parent_box;
    for dim = 1:3
        if bitget(c, dim)
            ch(1 + 2*dim) = md(dim); ch(2 + 2*dim) = hi(dim);
        else
            ch(1 + 2*dim) = lo(dim);  ch(2 + 2*dim) = md(dim);
        end
    end
    child_boxes(c+1, :) = ch;
end
end

function B = face_boxes(n)
% Floating-point chart endpoints for the 8 standard face charts of S^3.
% They are interpreted as INTLAB intervals in qn_single_box_certificate.m.
edges = linspace(-1, 1, n + 1);
B = zeros(8*n^3, 8);
r = 0;
for axisdim = 1:4
    for sgn = [1 -1]
        for i = 1:n
            for j = 1:n
                for k = 1:n
                    r = r + 1;
                    B(r, :) = [axisdim, sgn, edges(i), edges(i+1), ...
                               edges(j), edges(j+1), edges(k), edges(k+1)];
                end
            end
        end
    end
end
end
