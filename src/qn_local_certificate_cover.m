function qn_local_certificate_cover(worker_id, worker_count, face_subdivisions, outdir)
% Process one worker slice of the local single-box-certificate cover.
% This is the covering/subdivision loop around alg:single-box-certificate.
%
% PDF notation:
%   p = t e, e ∈ S³, 0 ≤ t ≤ ρ#                         Sec. 1
%   B = (E ∩ S³) × T                                    Eq. (42)
%   S(t,e), S_B                                         Eqs. (13), (67)
%   Single-box certificate returns VERIFIED/UNDECIDED
%
% Evaluation:
%   The face boxes E_j are ordinary floating-point endpoints, immediately
%   converted to INTLAB intervals inside qn_single_box_certificate.m. Each box is
%   certified by interval arithmetic; boxes failing only S(t,e)>0 are
%   recursively subdivided and re-tested.
%
% Writes per-box result lines to outdir/res_<procid>.csv and a DONE marker.
% Every real-valued CSV field is written with %.17e so the decimal strings
% retain round-trip double precision for the INTLAB interval endpoints.

if nargin < 4, outdir = 'results'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
PI = intval('pi');
rho = 3232/(27*PI^6);                         % rho# in Sec. 1
tmax = intval(sup(rho))*intval('1.000000000001');
radial_grid = tmax * (intval([0 2 4 6 8 10 12 14 15 16])/intval(16)); % T intervals, refined near ρ#
fid = fopen(sprintf('%s/res_%03d.csv', outdir, worker_id), 'w');
fprintf(fid, ['box_id,ok,inf_S,max_subdivision_depth,elapsed_seconds,' ...
    'lambda1_lower,lambda2_upper,mass_eigenvalue_lower,lambda3_lower,lambda3_gap_lower\n']);
face_box_table = face_boxes(face_subdivisions); % E boxes for B in Eq. (42)
total_boxes = size(face_box_table, 1);
tstart = tic;
failed_boxes = 0; worst_lower_S = inf;
for box_id = worker_id:worker_count:total_boxes
    face_box = face_box_table(box_id, :);       % one E; radial_grid supplies T
    proof = certify_with_subdivision(face_box, 0, radial_grid, [], [], []);
    worst_lower_S = min(worst_lower_S, proof.lower_S);
    if ~proof.verified, failed_boxes = failed_boxes + 1; end
    fprintf(fid, '%d,%d,%.17e,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n', ...
        box_id, proof.verified, proof.lower_S, proof.depth, toc(tstart), ...
        proof.lambda1_lower, proof.lambda2_upper,proof.mass_lower, ...
        proof.lambda3_lower,proof.lambda3_gap_lower);
    if mod(box_id, 40*worker_count) < worker_count
        fprintf('proc %d: box %d/%d worstS=%.3f fails=%d (%.0fs)\n', ...
                worker_id, box_id, total_boxes, worst_lower_S, failed_boxes, toc(tstart));
    end
end
fclose(fid);
fid2 = fopen(sprintf('%s/done_%03d.txt', outdir, worker_id), 'w');
fprintf(fid2, 'worstS=%.17e nfail=%d time=%.17e\n', worst_lower_S, failed_boxes, toc(tstart));
fclose(fid2);
end

function proof = certify_with_subdivision(face_box, depth, radial_grid, remainder_bounds, t_levels, parent_data)
maxdepth = 6;
if isempty(parent_data)
    % Full execution of alg:single-box-certificate.
    [res, child_remainder_bounds] = qn_single_box_certificate(face_box, [], radial_grid, remainder_bounds, t_levels);
    proof = summarize_results(res, depth);
    if proof.verified || depth >= maxdepth, return; end
    failed = find(~[res.ok]);
    sign_only = arrayfun(@(k) strcmp(res(k).reason,'S') && ...
        ~isempty(res(k).S_core_hessian), failed);
    sign_levels = failed(sign_only);
    child_parent_data = struct('REM',child_remainder_bounds, ...
        'S_padding',[res.S_padding], ...
        'veigs_verified',all([res(sign_levels).veigs_verified]));
    child_parent_data.S_core_hessian = {res.S_core_hessian};
    tasks = cell(0,3);
    if any(sign_only)
        tasks(end+1,:) = {sign_levels,child_remainder_bounds,child_parent_data};
    end
    if any(~sign_only)
        tasks(end+1,:) = {failed(~sign_only),child_remainder_bounds,[]};
    end
else
    % Cheap S-only re-test of the final sign condition S_B, Eq. (67).
    assert(isfield(parent_data,'veigs_verified') && parent_data.veigs_verified, ...
        'dq2:VEIGSInheritance', ...
        'A sign-only child may reuse only a parent box already verified by veigs.');
    res = dq2_retest_sign_box(face_box, radial_grid, t_levels, parent_data);
    proof = summarize_results(res, depth);
    if proof.verified || depth >= maxdepth, return; end
    if mod(depth,2)==1, next_parent=[]; else, next_parent=parent_data; end
    tasks = {find(~[res.ok]),parent_data.REM,next_parent};
end
proof.verified = true;
for child_box = subdivide_face_box(face_box)'
    for k=1:size(tasks,1)
        child = certify_with_subdivision(child_box',depth+1,radial_grid, ...
            tasks{k,2},tasks{k,1},tasks{k,3});
        proof = merge_proofs(proof,child);
    end
    if ~proof.verified && depth >= maxdepth-1, return; end
end
end

function proof = summarize_results(res,depth)
ok=[res.ok];
lower=min([res(ok).infS]);
if isempty(lower), lower=inf; end
proof=struct('verified',all(ok),'lower_S',lower,'depth',depth, ...
    'lambda1_lower',inf,'lambda2_upper',-inf,'mass_lower',inf, ...
    'lambda3_lower',inf,'lambda3_gap_lower',inf);
if ~isfield(res,'lam1'), return; end
for k=1:numel(res)
    if ~isempty(res(k).lam1)
        proof.lambda1_lower=min(proof.lambda1_lower,res(k).lam1(1));
        proof.lambda2_upper=max(proof.lambda2_upper,res(k).lam2(2));
    end
    if isfield(res,'mass_lower') && ~isempty(res(k).mass_lower)
        proof.mass_lower=min(proof.mass_lower,res(k).mass_lower);
        proof.lambda3_lower=min(proof.lambda3_lower,res(k).lambda3_lower);
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
