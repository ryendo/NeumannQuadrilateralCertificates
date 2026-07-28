function [verdict,info,gap,split_dim] = qn_certify_box(center,half_widths)
% Certified per-box index-1 eigenvalue bound using veigs.
% Floating-point work chooses only the subdivision coordinate.

C=qn_global_constants(); center=center(:); half_widths=half_widths(:);
[Kc,Mc]=qn_km_float(center);
[V,D]=eig((Kc+Kc')/2,(Mc+Mc')/2);
[evals,idx]=sort(real(diag(D))); V=real(V(:,idx));
gap=evals(2)-evals(1); v=V(:,1);

% Slack-driven coordinate, with longest-side fallback. This is non-certified.
slack=zeros(4,1);
for k=1:4
    h=1e-6*max(1,abs(center(k))); pp=center; pm=center;
    pp(k)=pp(k)+h; pm(k)=pm(k)-h;
    [Kp,Mp]=qn_km_float(pp); [Km,Mm]=qn_km_float(pm);
    dK=(Kp-Km)/(2*h); dM=(Mp-Mm)/(2*h);
    slack(k)=(abs(v'*dK*v)+abs(v'*dM*v))*half_widths(k);
end
[mx,split_dim]=max(slack);
[longest, longest_dim]=max(half_widths);
if ~(isfinite(mx) && mx>0) || half_widths(split_dim)<longest/2
    split_dim=longest_dim;
end

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
qh=intval(sup(area));

try
    [lambda1,veigs_info]=qn_veigs_smallest(K,M);
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
