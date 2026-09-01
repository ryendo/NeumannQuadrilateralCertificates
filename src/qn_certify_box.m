function [verdict,info,split_dim] = qn_certify_box(center,half_widths)
% Certified box test in Appendix B, equations (43)--(45).
% Floating-point calculations choose only the subdivision coordinate.  The
% certified decision uses all five conditions in (44).
% Here q_box=[q]_B encloses q(p)=|Q_p|; q_bounds contains the corresponding
% lower and upper endpoints used in (44).

C=qn_global_constants(); center=center(:); half_widths=half_widths(:);
split_dim=choose_split_dimension(center,half_widths);

try
    [K,M,q_box]=qn_km_enclosure(center,half_widths);
catch ME
    if startsWith(ME.identifier,'qn:')
        if strcmp(ME.identifier,'qn:QNotPositive')
            reason='q_not_positive';
        else
            reason=ME.identifier;
        end
        verdict='fail'; info=struct('reason',reason,'detail',ME.message); return;
    end
    rethrow(ME);
end
q_bounds=[inf(q_box),sup(q_box)];
if any(~isfinite(q_bounds)) || q_bounds(1)<=0
    verdict='fail'; info=struct('reason','q_not_positive', ...
        'q_lower',q_bounds(1),'q_upper',q_bounds(2),'q',q_bounds); return;
end
if ~qn_interval_ldl_pd(M)
    verdict='fail'; info=struct('reason','mass_not_pd', ...
        'q_lower',q_bounds(1),'q_upper',q_bounds(2),'q',q_bounds, ...
        'mass_pd',false); return;
end
try
    [bounds,veigs_info]=qn_veigs_indices(K,M,1);
    lambda1=bounds(1);
catch ME
    if startsWith(ME.identifier,'qn:VEIGS')
        verdict='fail'; info=struct('reason',ME.identifier,'detail',ME.message, ...
            'q_lower',q_bounds(1),'q_upper',q_bounds(2),'q',q_bounds, ...
            'mass_pd',true); return;
    end
    rethrow(ME);
end
index_range=veigs_info.indices;
lambda1_bounds=[inf(lambda1),sup(lambda1)];
[accepted,reason,delta_lower]= ...
    qn_global_code_test(q_box,true,lambda1,index_range,C.pi2);
if accepted
    verdict='cert';
    info=struct('route','veigs','delta_lower',delta_lower, ...
        'margin',delta_lower,'q_lower',q_bounds(1),'q_upper',q_bounds(2), ...
        'q',q_bounds,'mass_pd',true,'lambda1',lambda1_bounds, ...
        'index_range',index_range,'veigs_commit',veigs_info.commit);
else
    verdict='fail';
    info=struct('reason',reason,'delta_lower',delta_lower, ...
        'q_lower',q_bounds(1),'q_upper',q_bounds(2),'q',q_bounds, ...
        'mass_pd',true,'lambda1',lambda1_bounds,'index_range',index_range);
end
end

function split_dim=choose_split_dimension(center,half_widths)
% Coordinate choice described before (47).  It is a heuristic; the fallback
% guarantees h_r >= (1/2) max_s h_s, which is the condition used in the proof.
[longest,longest_dim]=max(half_widths);
split_dim=longest_dim;
try
    [Kc,Mc]=qn_km_float(center);
    [V,D]=eig((Kc+Kc')/2,(Mc+Mc')/2);
    [~,order]=sort(real(diag(D)));
    v=real(V(:,order(1)));
    if any(~isfinite(v)), return; end

    change=zeros(4,1);
    for k=1:4
        h=1e-6*max(1,abs(center(k)));
        p_plus=center; p_minus=center;
        p_plus(k)=p_plus(k)+h; p_minus(k)=p_minus(k)-h;
        [K_plus,M_plus]=qn_km_float(p_plus);
        [K_minus,M_minus]=qn_km_float(p_minus);
        dK=(K_plus-K_minus)/(2*h);
        dM=(M_plus-M_minus)/(2*h);
        change(k)=(abs(v'*dK*v)+abs(v'*dM*v))*half_widths(k);
    end
    [largest_change,candidate]=max(change);
    if isfinite(largest_change) && largest_change>0 && ...
            half_widths(candidate)>=longest/2
        split_dim=candidate;
    end
catch
    % Failure of this non-certified calculation leaves the longest side.
end
end
