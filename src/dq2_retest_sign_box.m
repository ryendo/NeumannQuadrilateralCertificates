function results = dq2_retest_sign_box(child_box, radial_grid, t_levels, parent_data)
% Cheap re-test of S-only failures on a sub-box.
% This only rechecks the final sign condition inf S_B > 0 from
% the final steps of Algorithm 1; the ancestor already certified
% Schur/root conditions.
%
% PDF notation:
%   S(t,e) is defined in (24), and positivity is required in (27).
%   Its coefficients come from the Vieta data d_2, tilde_d_1, d_0 in
%   (33)-(34), enclosed by (38) and optionally refined by (41).
%
% Evaluation:
%   Reuses ancestor INTLAB interval data for all non-S checks.  On a child
%   box it evaluates only the centered-form enclosure of S(t,e), using the
%   inherited Hessian bound and inherited root-error padding.
%
% Valid because every non-S check of the ancestor box (windows, localization,
% error constants S_padding, box Hessian of S_core) holds a fortiori on any sub-box.
% child_box: child face box. radial_grid: full t-grid. t_levels: levels to re-test.
% parent_data: ancestor data: .REM (remainder bounds), .S_padding(k),
%     .S_core_hessian{k} (Hessian objects over the ancestor box).

FR = dq2_fixed_reference();
TC = dq2_load_taylor_coefficients();

axisdim = child_box(1); sgn = child_box(2);
xI = intval(zeros(3,1));
for d = 1:3
    xI(d) = infsup(child_box(1+2*d), child_box(2+2*d));
end
xc = intval(mid(xI));
% child center must lie in the ancestor box for the centered form; it does,
% since the child is a sub-box. dxA = child box minus child center.
dx = xI - xc;

xhp = hessianinit(xc);
EP = dq2_face_direction(xhp, axisdim, sgn);
EK = midrad(zeros(5), parent_data.REM.BK); EM = midrad(zeros(5), parent_data.REM.BM);
Eb = midrad(zeros(5,1), parent_data.REM.Bb);
[CKP, CMP, CbP] = dq2_evaluate_taylor_coefficients_vectorized(TC, EP(1), EP(2), EP(3), EP(4)); % symbolic coefficients -> center Hessian data
qP = EP(1)*EP(1) + EP(4)*EP(4);

m = size(radial_grid, 2) - 1;
results = repmat(struct('ok',true,'infS',inf,'reason','skip'), m, 1);
for k = t_levels(:)'
    t = dq2_interval_hull(radial_grid(k), radial_grid(k+1));
    [~, d2_P, ~, S_core_center_data] = qn_single_box_quantities(CKP, CMP, CbP, EK, EM, Eb, qP, t, FR); % centered [S]_B data from (24), (38)
    if isempty(S_core_center_data) || isempty(d2_P), results(k).ok = false; results(k).infS = -inf; results(k).reason = 'cinv'; continue; end
    % centered form: point value+gradient at child center, ancestor Hessian
    S_core_hessian = parent_data.S_core_hessian{k};
    center_value = S_core_center_data.x; center_gradient = S_core_center_data.dx;
    S_centered = center_value;
    for i = 1:3
        S_centered = S_centered + center_gradient(i)*dx(i);
    end
    quadratic_part = intval(0);
    H = S_core_hessian.hx;
    for i = 1:3
        for j = 1:3
            quadratic_part = quadratic_part + H(i,j)*dx(i)*dx(j);
        end
    end
    S_core_child = S_centered + quadratic_part/2;
    S_B_child = S_core_child + midrad(0, parent_data.S_padding(k)); % final sign test from (24), (27)
    results(k).infS = inf(S_B_child);
    results(k).ok = inf(S_B_child) > 0;
    if ~results(k).ok, results(k).reason = 'S'; end
end
end
