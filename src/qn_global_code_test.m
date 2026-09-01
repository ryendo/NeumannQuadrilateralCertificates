function [accepted,reason,delta_lower] = ...
    qn_global_code_test(q_box,mass_pd,lambda_box,index_range,pi2)
% The five accepted-box conditions displayed in equation (44).  In the
% paper's notation q_box=[q]_B and its endpoints are lower(q_B), upper(q_B);
% q denotes area, while A denotes the pulled-back form matrix of Section 4.3.

accepted=false; delta_lower=nan;
q_lower=inf(q_box); q_upper_endpoint=sup(q_box);
if ~isscalar(q_lower) || ~isscalar(q_upper_endpoint) || ...
        ~isfinite(q_lower) || ~isfinite(q_upper_endpoint) || q_lower<=0
    reason='q_not_positive';
    return
end
if ~isscalar(mass_pd) || mass_pd~=true
    reason='mass_not_pd';
    return
end
if ~isnumeric(index_range) || ~any(index_range(:)==1)
    reason='index_1_not_verified';
    return
end
lambda_lower=inf(lambda_box); lambda_upper_endpoint=sup(lambda_box);
if ~isscalar(lambda_lower) || ~isscalar(lambda_upper_endpoint) || ...
        ~isfinite(lambda_lower) || ~isfinite(lambda_upper_endpoint) || ...
        lambda_lower<=0
    reason='lambda1_not_positive';
    return
end

pi2_lower=inf(pi2);
if ~isscalar(pi2_lower) || ~isfinite(pi2_lower)
    reason='delta_not_positive';
    return
end
q_upper=intval(q_upper_endpoint);
lambda_upper=intval(lambda_upper_endpoint);
delta_lower=inf(intval(pi2_lower)-q_upper*lambda_upper);
if ~isfinite(delta_lower) || delta_lower<=0
    reason='delta_not_positive';
    return
end
accepted=true; reason='accepted';
end
