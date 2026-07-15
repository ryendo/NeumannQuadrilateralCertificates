function dq2_run_certificate(worker_id, worker_count, face_subdivisions, ~, ~, ~, outdir)
% Process one worker slice of the S^3 face-box covering.
% This is the outer covering/subdivision loop around Algorithm 1 in Sec. 4.1.
%
% PDF notation:
%   p = t e, e ∈ S³, 0 ≤ t ≤ ρ#                         Sec. 1
%   B = (E ∩ S³) × T                                    Eq. (42)
%   S(t,e), S_B                                         Eqs. (13), (67)
%   Algorithm 1 returns VERIFIED/UNDECIDED              Sec. 4.1
%
% Evaluation:
%   The face boxes E_j are ordinary floating-point endpoints, immediately
%   converted to INTLAB intervals inside dq2_algorithm1_box.m.  Each box is
%   certified by interval arithmetic; boxes failing only S(t,e)>0 are
%   recursively subdivided and re-tested.
%
% Writes per-box result lines to outdir/res_<procid>.csv and a DONE marker.
% Every real-valued CSV field is written with %.17e so the decimal strings
% retain round-trip double precision for the INTLAB interval endpoints.
% m is kept for interface compatibility; the t-grid refines the top levels.

if nargin < 7, outdir = 'results'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
PI = intval('pi');
rho = 3232/(27*PI^6);                         % rho# in Sec. 1
tmax = intval(sup(rho))*intval('1.000000000001');
radial_grid = tmax * (intval([0 2 4 6 8 10 12 14 15 16])/intval(16)); % T intervals, refined near ρ#
fid = fopen(sprintf('%s/res_%03d.csv', outdir, worker_id), 'w');
fprintf(fid, 'box_id,ok,inf_S,max_subdivision_depth,elapsed_seconds,lambda1_lower,lambda2_upper\n');
face_box_table = face_boxes(face_subdivisions); % E boxes for B in Eq. (42)
total_boxes = size(face_box_table, 1);
tstart = tic;
failed_boxes = 0; worst_lower_S = inf;
for box_id = worker_id:worker_count:total_boxes
    face_box = face_box_table(box_id, :);       % one E; radial_grid supplies T
    [is_verified, lower_S, max_depth, lambda1_lower, lambda2_upper] = ...
        certify_with_subdivision(face_box, 0, radial_grid, [], [], []);
    worst_lower_S = min(worst_lower_S, lower_S);
    if ~is_verified, failed_boxes = failed_boxes + 1; end
    fprintf(fid, '%d,%d,%.17e,%d,%.17e,%.17e,%.17e\n', ...
            box_id, is_verified, lower_S, max_depth, toc(tstart), lambda1_lower, lambda2_upper);
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

function [is_verified, lower_S, depth_used, lambda1_lower, lambda2_upper] = certify_with_subdivision(face_box, depth, radial_grid, remainder_bounds, t_levels, parent_data)
maxdepth = 6;
depth_used = depth;
lambda1_lower = inf; lambda2_upper = -inf;
if isempty(parent_data)
    % Full single-box certificate: Algorithm 1, Steps 1-14.
    [res, child_remainder_bounds] = dq2_algorithm1_box(face_box, [], radial_grid, remainder_bounds, t_levels);
    for k = 1:numel(res)
        if ~isempty(res(k).lam1)
            lambda1_lower = min(lambda1_lower, res(k).lam1(1));
            lambda2_upper = max(lambda2_upper, res(k).lam2(2));
        end
    end
    okv = [res.ok];
    lower_S = min([res(okv).infS]);
    if isempty(lower_S), lower_S = inf; end
    is_verified = all(okv);
    if is_verified || depth >= maxdepth, return; end
    sfail = []; wfail = [];
    for k = find(~okv)
        if strcmp(res(k).reason, 'S') && ~isempty(res(k).S_core_hessian)
            sfail(end+1) = k; %#ok<AGROW>
        else
            wfail(end+1) = k; %#ok<AGROW>
        end
    end
    child_parent_data = struct( ...
        'REM', child_remainder_bounds, ...
        'S_padding', [res.S_padding], ...
        'S_core_hessian', {{res.S_core_hessian}}); % interval data inherited by child boxes
    child_parent_data.S_core_hessian = {res.S_core_hessian};
    is_verified = true;
    for child_box = subdivide_face_box(face_box)'
        if ~isempty(sfail)
            [child_verified, child_lower_S, child_depth, child_lambda1_lower, child_lambda2_upper] = ...
                certify_with_subdivision(child_box', depth+1, radial_grid, child_remainder_bounds, sfail, child_parent_data);
            is_verified = is_verified && child_verified; lower_S = min(lower_S, child_lower_S); depth_used = max(depth_used, child_depth);
            lambda1_lower = min(lambda1_lower, child_lambda1_lower); lambda2_upper = max(lambda2_upper, child_lambda2_upper);
        end
        if ~isempty(wfail)
            [child_verified, child_lower_S, child_depth, child_lambda1_lower, child_lambda2_upper] = ...
                certify_with_subdivision(child_box', depth+1, radial_grid, child_remainder_bounds, wfail, []);
            is_verified = is_verified && child_verified; lower_S = min(lower_S, child_lower_S); depth_used = max(depth_used, child_depth);
            lambda1_lower = min(lambda1_lower, child_lambda1_lower); lambda2_upper = max(lambda2_upper, child_lambda2_upper);
        end
        if ~is_verified && depth >= maxdepth - 1, return; end
    end
else
    % Cheap S-only re-test of the final sign condition S_B, Eq. (67).
    res = dq2_retest_sign_box(face_box, radial_grid, t_levels, parent_data);
    okv = [res.ok];
    lower_S = min([res(okv).infS]);
    if isempty(lower_S), lower_S = inf; end
    is_verified = all(okv);
    if is_verified || depth >= maxdepth, return; end
    failedlv = find(~okv);
    is_verified = true;
    for child_box = subdivide_face_box(face_box)'
        if mod(depth, 2) == 1
            % refresh with a full evaluation every other depth
            [child_verified, child_lower_S, child_depth, child_lambda1_lower, child_lambda2_upper] = ...
                certify_with_subdivision(child_box', depth+1, radial_grid, parent_data.REM, failedlv, []);
        else
            [child_verified, child_lower_S, child_depth, child_lambda1_lower, child_lambda2_upper] = ...
                certify_with_subdivision(child_box', depth+1, radial_grid, parent_data.REM, failedlv, parent_data);
        end
        is_verified = is_verified && child_verified; lower_S = min(lower_S, child_lower_S); depth_used = max(depth_used, child_depth);
        lambda1_lower = min(lambda1_lower, child_lambda1_lower); lambda2_upper = max(lambda2_upper, child_lambda2_upper);
        if ~is_verified && depth >= maxdepth - 1, return; end
    end
end
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
% They are interpreted as INTLAB intervals in dq2_algorithm1_box.m.
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
