function [verdict,info,split_dim] = qn_certify_box(center,half_widths)
% Certified per-box index-1 eigenvalue bound using veigs.
% An undecided box is bisected along a longest side.

C=qn_global_constants(); center=center(:); half_widths=half_widths(:);
[~,split_dim]=max(half_widths);

try
    [K,M,area]=qn_km_enclosure(center,half_widths);
catch ME
    if startsWith(ME.identifier,'qn:')
        verdict='fail'; info=struct('reason',ME.identifier); return;
    end
    rethrow(ME);
end
if ~qn_interval_ldl_pd(M)
    verdict='fail'; info=struct('reason','mass_not_pd'); return;
end
qh=intval(sup(area)); % paper's qbar_B = sup_{p in B} q(p) = sup_{p in B}|Q_p|

try
    [bounds,veigs_info]=qn_veigs_indices(K,M,1);
    lambda1=bounds(1);
catch ME
    if startsWith(ME.identifier,'qn:VEIGS')
        verdict='fail'; info=struct('reason',ME.identifier,'detail',ME.message); return;
    end
    rethrow(ME);
end
ub=intval(sup(lambda1));
if sup(qh*ub)<inf(C.pi2)
    verdict='cert';
    info=struct('route','veigs','margin',inf(C.pi2-qh*ub), ...
        'lambda1',[inf(lambda1),sup(lambda1)], ...
        'index_range',veigs_info.indices,'veigs_commit',veigs_info.commit);
else
    verdict='fail';
    info=struct('reason','veigs_bound','lambda1',[inf(lambda1),sup(lambda1)], ...
        'index_range',veigs_info.indices);
end
end
