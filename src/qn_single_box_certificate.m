function results = qn_single_box_certificate(face_box, radial_grid, t_levels)
% Certified INTLAB second-order difference-quotient test on one face box.
% Implements Algorithm 1, "Certificate on B," in Appendix A.4.
%
% Current paper formulas:
%   λᵢ(te) = π² + tνᵢ(t,e), and L is the scaled cluster sum.           (22)
%   ν₁(t,e)ν₂(t,e) has the continuous extension displayed before (23).
%   S(t,e) is the expression in ν₁ν₂ and L given in (23).
%   C_t(ν) is the integral difference-quotient pencil.                (27)
%   Φ_t(ν,e) and F_t(ν) define the 2×2 Schur complement.              (28)
%   det F_t(ν)=0 is equivalent to the full generalized eigenproblem.  (30)
%   F̂_t and Ψ_t are defined before (31); their comparison identity is (33).
%   d̃₁ is defined in (31), and det F̂_t is written in (32).
%   L and ν₁ν₂ are enclosed by (36), with the refinement in (38)--(39).
%
% Code variable map:
%   face_box                  -> chart box from χ_{r,σ} in (40)
%   t                         -> radial interval T = [t₋,t₊]
%   I_nu, I_nu_tight          -> interval [ν₋,ν₊] chosen in Algorithm 1
%   d2_B, d1_over_t_B, d0_B   -> d₂, d̃₁, d₀ in (31)-(32)
%   S_core_B, S_padding       -> central part and error radius of [S]_B
%                                obtained from (23) and (36)
%   S_B                       -> final [S]_B enclosure in Algorithm 1
%   root_error_radius         -> tη/β from the proof of (36)
%   mass_lower                -> Gershgorin lower bound for lambda_min(M(te))
%   lam1, lam3                -> verified index-1/index-3 enclosures
%
% Evaluation:
%   All rigorous quantities are INTLAB intervals.  Taylor coefficients of
%   K, M, and b are precomputed symbolically and loaded by
%   dq2_load_taylor_coefficients.m; they are then evaluated either as INTLAB
%   intervals (whole box) or as INTLAB Hessian objects (centered form).
%
% face_box = [axisdim sgn lo1 hi1 lo2 hi2 lo3 hi3].
% t_levels is an optional subset of the intervals in radial_grid.

[FR,K0,M0] = dq2_fixed_reference();
PI2 = FR.PI2;
TC = dq2_load_taylor_coefficients();

axisdim = face_box(1); sgn = face_box(2);
xI = intval(zeros(3,1));
for d = 1:3
    xI(d) = infsup(face_box(1+2*d), face_box(2+2*d));
end
xc = intval(mid(xI));

if nargin < 2 || isempty(radial_grid)
    error('dq2:RadialGrid','The certified radial grid is required.');
end
radial_edges = radial_grid;
m = size(radial_edges, 2) - 1;

% Direction enclosures for B = (E intersect S^3) x T, Eq. (40).
EI = dq2_face_direction(xI, axisdim, sgn);   % e-box E_j, evaluated as INTLAB intervals
xhp = hessianinit(xc);
EP = dq2_face_direction(xhp, axisdim, sgn);  % e at center, evaluated as INTLAB Hessian variables
xhb = hessianinit(xI);
EG = dq2_face_direction(xhb, axisdim, sgn);  % e over E_j, Hessian enclosure for centered form

% Taylor-quotient remainder bounds implement the integral in C_t from (27),
% avoiding division by an interval containing t = 0.
[BKr, BMr, Bbr] = dq2_bound_taylor_remainder_vectorized(EI, radial_edges(end), 2, 10);
EK = midrad(zeros(5), BKr); EM = midrad(zeros(5), BMr); Eb = midrad(zeros(5,1), Bbr); % interval R_K, R_M, R_b
if nargin < 3 || isempty(t_levels)
    t_levels = 1:m;
end

% Symbolic Taylor coefficient table evaluated in three rigorous modes:
[CKI, CMI, CbI] = dq2_evaluate_taylor_coefficients_vectorized(TC, EI(1), EI(2), EI(3), EI(4)); % whole-box INTLAB intervals
[CKP, CMP, CbP] = dq2_evaluate_taylor_coefficients_vectorized(TC, EP(1), EP(2), EP(3), EP(4)); % center value/gradient
[CKG, CMG, CbG] = dq2_evaluate_taylor_coefficients_vectorized(TC, EG(1), EG(2), EG(3), EG(4)); % Hessian enclosure over box

ad2I = sqr(EI(1)) + sqr(EI(4));               % e_a^2+e_d^2 in (23), not the area q(p)
ad2P = EP(1)*EP(1) + EP(4)*EP(4);
ad2G = EG(1)*EG(1) + EG(4)*EG(4);

results = repmat(struct('ok',true,'infS',inf,'reason','skip','detail','','lam1',[],'lam2',[], ...
                        'lam3',[],'mass_lower',-inf,'lambda3_gap_lower',-inf, ...
                        'exact_refinement_used',false, ...
                        'S_before_exact_refinement',-inf, ...
                        'exact_refinement_gain',0,'veigs_indices',[]), m, 1);
for k = t_levels(:)'
    results(k).ok = false; results(k).infS = -inf; results(k).reason = '';
    t = dq2_interval_hull(radial_edges(k), radial_edges(k+1));
    % -- Inverse-free precheck: raw Taylor blocks, M, and B0 at nu=0 -----
    blocksI = qn_single_box_blocks(CKI, CMI, CbI, EK, EM, Eb, ad2I, t, FR);
    Mte = M0 + t*blocksI.DtM;                 % M(te), positive-definite condition in (34)
    Kte = K0 + t*blocksI.DtK;                 % K(te), from the same Taylor quotient
    mass_lower = gershgorin_lower_bound(Mte);
    if ~(isfinite(mass_lower) && mass_lower > 0)
        results(k).reason = 'mass_not_pd'; continue;
    end
    if ~posdef3(blocksI.B0)                   % D_perp+t*C_{t,perp,perp}(0)
        results(k).reason = 'affine_perp_not_pd'; continue;
    end
    % -- Affine Schur data: B0 is inverted only after the PD precheck -----
    [d1_over_t_I, d2_I, d0_I, S_core_I, schur] = qn_single_box_quantities(blocksI, ad2I, t, FR); % whole-box interval evaluation
    if isempty(schur), results(k).reason = 'affine_perp_inverse'; continue; end
    % ---------- centered chains (hessian at center / over box) ----------
    blocksP = qn_single_box_blocks(CKP, CMP, CbP, EK, EM, Eb, ad2P, t, FR);
    blocksG = qn_single_box_blocks(CKG, CMG, CbG, EK, EM, Eb, ad2G, t, FR);
    [d1_over_t_P, d2_P, d0_P, S_core_P] = qn_single_box_quantities(blocksP, ad2P, t, FR); % center Hessian object
    [d1_over_t_G, d2_G, d0_G, S_core_G] = qn_single_box_quantities(blocksG, ad2G, t, FR); % box Hessian enclosure
    if isempty(d1_over_t_P) || isempty(d1_over_t_G), results(k).reason = 'centered_inverse'; continue; end
    dx = xI - xc;
    d1_over_t_B = intersect(centered_hessian_interval(d1_over_t_P, d1_over_t_G, dx), d1_over_t_I); % d̃₁, (31)-(32)
    d2_B = intersect(centered_hessian_interval(d2_P, d2_G, dx), d2_I);                             % d₂, (32)
    d0_B = intersect(centered_hessian_interval(d0_P, d0_G, dx), d0_I);                             % d₀, (32)
    S_core_B = intersect(centered_hessian_interval(S_core_P, S_core_G, dx), S_core_I);             % central [S]_B part from (23), (36)
    if isnan(inf(d1_over_t_B)) || isnan(inf(d2_B)) || isnan(inf(d0_B)) || isnan(inf(S_core_B))
        results(k).reason = 'centered_intersection'; continue;
    end
    if inf(d2_B) <= 0, results(k).reason = 'affine_root_leading_coefficient'; continue; end
    % Intervals for the roots nuhat_1,nuhat_2 of (32).
    d1_B = t*d1_over_t_B;
    DiscC = max(sqr(d1_B) - 4*d2_B*d0_B, intval(0));
    rt = sqrt(DiscC);
    nuhat1 = (d1_B - rt)/(2*d2_B);
    nuhat2 = (d1_B + rt)/(2*d2_B);
    beta = inf(lammin2(schur.minus_dFhat_dnu)); % beta in (35)
    if beta <= 0, results(k).reason = 'affine_slope_not_pd'; continue; end
    pad = upper_initial_pad(schur.al, schur.be);
    ok = true; reason = ''; eta = inf; root_error_radius = inf; I_nu = [];
    nw = 6;   % subdivision of [nu] used to certify (34)--(35)
    L1 = -schur.Nm + schur.C0*schur.B0i*schur.MW;
    nrmL1 = spnorm(L1);
    window_closed = false;
    for it = 1:3
        I_nu = padded_hull(nuhat1, nuhat2, pad);
        [kappa2, ok, reason] = schur_residual_norm_bound(schur, t, I_nu, nw, nrmL1); % interval bound for ||F_t - Fhat_t||/t
        if ~ok, break; end
        eta = upper_eta(I_nu, t, kappa2);       % eta in (35)
        root_error_radius = upper_root_error_radius(t, eta, beta); % t*eta/beta in the proof of (36)
        if upper_plus_margin(root_error_radius) <= pad
            window_closed = true;
            break
        end
        pad = upper_scaled_margin(root_error_radius);
    end
    if ~ok, results(k).reason = reason; continue; end
    if ~window_closed, results(k).reason = 'nu_interval_not_closed'; continue; end
    % The roots lie in hull(nuhat1,nuhat2) +/- root_error_radius; re-evaluate the
    % comparison bound on this smaller interval.
    for bs = 1:2
        I_nu_tight = intersect(padded_hull(nuhat1, nuhat2, root_error_radius), I_nu);
        [kapr, ok, reason] = schur_residual_norm_bound(schur, t, I_nu_tight, nw, nrmL1);
        if ~ok, break; end
        eta_n = upper_eta(I_nu_tight, t, kapr);
        radius_n = upper_root_error_radius(t, eta_n, beta);
        if radius_n >= root_error_radius, break; end
        eta = eta_n; root_error_radius = radius_n;
    end
    if ~ok, results(k).reason = reason; continue; end
    ok = true;
    Xfull = schur.X0 + t*schur.Xt;            % E_0 block of the exact Schur pencil, interval
    wedges = interval_edges(I_nu, nw);
    for w = 1:nw
        Wj = dq2_interval_hull(wedges(w), wedges(w+1));
        BWj = schur.B0 - (t*Wj)*schur.MW;     % D_perp+t*C_{t,perp,perp}(nu), (34)
        if ~posdef3(BWj), ok = false; reason = 'perp_block_not_pd'; break; end
        [BWji, okinv] = inv3(BWj);
        if ~okinv, ok = false; reason = 'perp_block_inverse'; break; end
        CWj = schur.C0 - (t*Wj)*schur.Nm;
        dZ = -(schur.Nm*BWji*CWj' + CWj*BWji*schur.Nm') + CWj*BWji*schur.MW*BWji*CWj';
        mono = schur.Y + sqr(t)*symm(dZ);     % -partial_nu F_t(nu), (34)
        if inf(lammin2(symm(mono))) <= 0, ok = false; reason = 'schur_derivative'; break; end
    end
    if ~ok, results(k).reason = reason; continue; end
    [F_at_nu_minus,okminus] = schur_matrix_at_nu(Xfull, schur.Y, schur.C0, schur.B0, schur.Nm, schur.MW, t, intval(inf(I_nu))); % F_t(nu_-), (34)
    [F_at_nu_plus,okplus] = schur_matrix_at_nu(Xfull, schur.Y, schur.C0, schur.B0, schur.Nm, schur.MW, t, intval(sup(I_nu))); % F_t(nu_+), (34)
    if ~(okminus && okplus), results(k).reason = 'perp_block_endpoint'; continue; end
    if ~(inf(F_at_nu_minus(1,1)) > 0 && inf(det2(F_at_nu_minus)) > 0), results(k).reason = 'left_schur_not_pd'; continue; end
    if ~(sup(F_at_nu_plus(1,1)) < 0 && inf(det2(F_at_nu_plus)) > 0), results(k).reason = 'right_schur_not_nd'; continue; end
    sbar = upper_abs(I_nu);
    L_padding = sup(intval(2)*intval(eta)/intval(beta)); % L enclosure radius in (36)
    nu_product_padding = sup(intval(2)*intval(sbar)*intval(root_error_radius) + ...
        sqr(intval(root_error_radius))); % radius for ν₁ν₂ in (36)
    S_padding = sup(intval(L_padding)*PI2*(intval(1) + intval(2)*sqr(t)) + ...
               intval(2)*intval(nu_product_padding)*(intval(1) + sqr(t)));
    S_B = S_core_B + midrad(0, S_padding);    % [S]_B from (23) and (36)
    nu1_interval = nuhat1 + midrad(0, root_error_radius); % [ν₁] in Algorithm 1
    nu2_interval = nuhat2 + midrad(0, root_error_radius); % [ν₂] in Algorithm 1
    % Retain the enclosure from (36) before intersecting with (38)--(39).
    S_before_exact_refinement = S_B;
    exact_refinement_used = false;
    % ---- exact refinement (38)--(39), used only for disjoint intervals ----
    gap = nu2_interval - nu1_interval;
    if inf(gap) > 0
        Fhat_at_zero = schur.X0 + t*schur.G;
        root_shift_over_t = intval(zeros(2,1)); exact_refinement_ok = true;
        nu_intervals = {nu1_interval, nu2_interval};
        for kk = 1:2
            nu_k = nu_intervals{kk};
            B_at_nu_k = schur.B0 - (t*nu_k)*schur.MW;
            if ~posdef3(B_at_nu_k), exact_refinement_ok = false; break; end
            [B_at_nu_k_inv, okinv] = inv3(B_at_nu_k);
            if ~okinv, exact_refinement_ok = false; break; end
            L2_at_nu_k = -B_at_nu_k_inv*schur.Nm' + schur.B0i*schur.MW*B_at_nu_k_inv*schur.C0';
            L1L2_at_nu_k = L1*L2_at_nu_k;
            L1L2_at_nu_k = (L1L2_at_nu_k + L1L2_at_nu_k')/2;
            residual_over_t = sqr(t)*sqr(nu_k)*L1L2_at_nu_k;       % E(ν)/t with E = t³ν²L₁L₂
            affine_at_nu_k = Fhat_at_zero - ...
                nu_k*schur.minus_dFhat_dnu; % Fhat_t(nu_k), (38)
            adj_affine = [affine_at_nu_k(2,2), -affine_at_nu_k(1,2); -affine_at_nu_k(2,1), affine_at_nu_k(1,1)];
            trace_adj_residual = adj_affine(1,1)*residual_over_t(1,1) + adj_affine(1,2)*residual_over_t(2,1) + ...
                                 adj_affine(2,1)*residual_over_t(1,2) + adj_affine(2,2)*residual_over_t(2,2);
            det_residual_over_t = t*sqr(sqr(t))*sqr(sqr(nu_k))*det2(L1L2_at_nu_k);   % det(E)/t = t⁵ν⁴det(L₁L₂)
            if kk == 1, root_gap = nu1_interval - nu2_interval; else, root_gap = nu2_interval - nu1_interval; end
            root_shift_over_t(kk) = (trace_adj_residual - det_residual_over_t)/(d2_B*root_gap);
        end
        if exact_refinement_ok
            exact_refinement_used = true;
            [L_correction,nu_product_correction] = ...
                qn_exact_root_refinement_corrections( ...
                nuhat1,nuhat2,root_shift_over_t(1),root_shift_over_t(2),t);
            S_correction = -PI2*L_correction - 2*nu_product_correction + ...
                2*ad2I*sqr(t)*(PI2*L_correction + nu_product_correction); % (39) substituted in (23)
            S_B_from_exact_refinement = S_core_B + S_correction;
            S_B = intersect(S_B, S_B_from_exact_refinement);
            nu1_interval = intersect(nu1_interval, nuhat1 + t*root_shift_over_t(1));
            nu2_interval = intersect(nu2_interval, nuhat2 + t*root_shift_over_t(2));
        end
    end
    if isnan(inf(S_B)) || isnan(sup(S_B)) || ...
            isnan(inf(nu1_interval)) || isnan(sup(nu1_interval)) || ...
            isnan(inf(nu2_interval)) || isnan(sup(nu2_interval))
        results(k).reason = 'exact_refinement_intersection'; continue;
    end
    lam1 = PI2 + t*nu1_interval;              % refined λ₁(te), (22)
    lam2 = PI2 + t*nu2_interval;              % refined λ₂(te), (22)
    % Independently target indices 1 and 3 of the complete 5x5 interval
    % pencil with verified veigs/veig calls. The second-order
    % difference-quotient analysis remains necessary for the double cluster
    % and S(t,e); acceptance additionally requires
    % the verified separation lambda_3 > 16.
    try
        [veigs_bounds,veigs_info]=qn_veigs_indices(Kte,Mte,[1 3]);
        veigs_lam1=veigs_bounds(1);
        veigs_lam3=veigs_bounds(2);
    catch ME
        if startsWith(ME.identifier,'qn:VEIGS')
            results(k).reason=ME.identifier; results(k).detail=ME.message; continue;
        end
        rethrow(ME);
    end
    results(k).mass_lower=mass_lower;
    results(k).lam3=[inf(veigs_lam3),sup(veigs_lam3)];
    results(k).veigs_indices=veigs_info.indices;
    lam1=intersect(lam1,veigs_lam1);
    if isnan(inf(lam1)) || isnan(sup(lam1))
        results(k).reason='veigs_intersection'; continue;
    end
    if inf(veigs_lam3) <= 16
        results(k).reason='veigs_lambda3'; continue;
    end
    results(k).infS = inf(S_B);
    results(k).lam1 = [inf(lam1), sup(lam1)];
    results(k).lam2 = [inf(lam2), sup(lam2)];
    results(k).lambda3_gap_lower = inf(veigs_lam3 - intval('1.5')*PI2);
    results(k).exact_refinement_used = exact_refinement_used;
    results(k).S_before_exact_refinement = inf(S_before_exact_refinement);
    results(k).exact_refinement_gain = ...
        inf(S_B) - inf(S_before_exact_refinement);
    if inf(S_B) <= 0
        results(k).reason = 'S_not_positive';
    elseif inf(lam1) <= 0
        % This is not an S-only failure: subdivision must rerun the complete
        % certificate rather than inherit the parent spectral enclosure.
        results(k).reason = 'lambda1';
    elseif mass_lower > 0 && inf(veigs_lam3) > 16
        results(k).ok = true;
    else
        results(k).reason = 'spectral';
    end
end
end

% ---------- helpers ----------
function c = centered_hessian_interval(cp, cb, dx)
% Centered-form interval evaluation: value and gradient at the center point
% plus Hessian enclosure over the box. Used for d_2, tilde_d_1, d_0, and S.
v = cp.x;                    % intval value at center
g = cp.dx;                   % 1x3 gradient at center
H = cb.hx;                   % 3x3 Hessian enclosure over the box
c = v;
for i = 1:3
    c = c + g(i)*dx(i);
end
qd = intval(0);
for i = 1:3
    for j = 1:3
        qd = qd + H(i,j)*dx(i)*dx(j);
    end
end
c = c + qd/2;
end

function S = symm(A)
S = (A + A')/2;
end

function lower = gershgorin_lower_bound(A)
% Certified lower bound for lambda_min(A) from symmetric Gershgorin discs.
% mag(...) is converted back to an interval before summation so every row
% radius is rounded upward.
A = hull(intval(A),intval(A)');
lower = inf;
for i = 1:size(A,1)
    jj = [1:i-1, i+1:size(A,2)];
    off = sum(intval(mag(A(i,jj))));
    lower = min(lower,inf(A(i,i)-off));
end
end

function d = det2(A)
d = A(1,1)*A(2,2) - A(1,2)*A(2,1);
end

function [Ai, ok] = inv3(A)
D = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) ...
  - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) ...
  + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1));
if isintval(D)
    ok = (inf(D) > 0) || (sup(D) < 0);
else
    Dx = D.x;
    ok = (inf(Dx) > 0) || (sup(Dx) < 0);
end
if ~ok, Ai = A; return; end
Ai = [ A(2,2)*A(3,3)-A(2,3)*A(3,2), A(1,3)*A(3,2)-A(1,2)*A(3,3), A(1,2)*A(2,3)-A(1,3)*A(2,2);
       A(2,3)*A(3,1)-A(2,1)*A(3,3), A(1,1)*A(3,3)-A(1,3)*A(3,1), A(1,3)*A(2,1)-A(1,1)*A(2,3);
       A(2,1)*A(3,2)-A(2,2)*A(3,1), A(1,2)*A(3,1)-A(1,1)*A(3,2), A(1,1)*A(2,2)-A(1,2)*A(2,1) ];
Ai = Ai/D;
end

function l = lammin2(A)
tr = A(1,1) + A(2,2);
disc = max(sqr(A(1,1)-A(2,2)) + 4*sqr(A(1,2)), intval(0));
l = (tr - sqrt(disc))/2;
end

function ok = posdef3(A)
m2 = A(1,1)*A(2,2) - A(1,2)*A(2,1);
d3 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) ...
   - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) ...
   + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1));
ok = inf(A(1,1)) > 0 && inf(m2) > 0 && inf(d3) > 0;
end

function edges = interval_edges(X, n)
% Interval endpoints for a subdivision of X; avoids double linspace.
lo = intval(inf(X)); hi = intval(sup(X));
edges = intval(zeros(1, n+1));
for j = 0:n
    edges(j+1) = lo + (hi - lo)*(intval(j)/intval(n));
end
edges(1) = lo;
edges(end) = hi;
end

function r = upper_abs(X)
% Upward-rounded scalar bound for max(abs(X)).
r = sup(abs(intval(X)));
end

function pad = upper_initial_pad(al, be)
% Adaptive root-window padding, computed as an upward-rounded interval bound.
pad = sup(intval(1) + intval(2)*(intval(rad(al)) + intval(rad(be))));
end

function W = padded_hull(nuhat1, nuhat2, pad)
% [inf(nuhat1)-pad, sup(nuhat2)+pad], with INTLAB endpoint arithmetic.
lo = inf(intval(inf(nuhat1)) - intval(pad));
hi = sup(intval(sup(nuhat2)) + intval(pad));
W = infsup(lo, hi);
end

function eta = upper_eta(W, t, kappa2)
% Upward-rounded eta bound for (35).
sbar = upper_abs(W);
eta = sup(sqr(intval(sbar)) * sqr(intval(t)) * intval(kappa2));
end

function radius = upper_root_error_radius(t, eta, beta)
% Upward-rounded root-error radius t*eta/beta from the proof of (36).
radius = sup(intval(t) * intval(eta) / intval(beta));
end

function y = upper_plus_margin(x)
y = sup(intval(x) + intval('1e-6'));
end

function y = upper_scaled_margin(x)
y = sup(intval('1.5') * (intval(x) + intval('1e-6')));
end

function n = spnorm(A)
% upper bound on ||A||_2 via lammax of the 2x2 Gram matrix
A = intval(A);
if size(A,1) <= size(A,2)
    Gm = A*A';
else
    Gm = A'*A;
end
tr = Gm(1,1) + Gm(2,2);
disc = max(sqr(Gm(1,1)-Gm(2,2)) + 4*sqr(Gm(1,2)), intval(0));
lmax = (tr + sqrt(disc))/2;
n = sup(sqrt(max(lmax, intval(0))));
end

function [kappa2, ok, reason] = schur_residual_norm_bound(schur, t, I_nu, nw, nrmL1)
% kappa2 = ||L1|| * sup over [nu] of ||L2(nu)||, using its subdivision.
wedges = interval_edges(I_nu, nw);
kappa2 = 0; ok = true; reason = '';
for w = 1:nw
    Wj = dq2_interval_hull(wedges(w), wedges(w+1));
    BWj = schur.B0 - (t*Wj)*schur.MW;
    if ~posdef3(BWj), ok = false; reason = 'perp_block_not_pd'; return; end
    [BWji, okinv] = inv3(BWj);
    if ~okinv, ok = false; reason = 'perp_block_inverse'; return; end
    L2j = -BWji*schur.Nm' + schur.B0i*schur.MW*BWji*schur.C0';
    kappa2 = max(kappa2, sup(intval(nrmL1)*intval(spnorm(L2j))));
end
end

function [F,ok] = schur_matrix_at_nu(Xfull, Y, C0, B0, Nm, MW, t, nu)
Cm = C0 - (t*nu)*Nm;
Bm = B0 - (t*nu)*MW;
ok = posdef3(Bm);
if ~ok
    F = intval(zeros(2));
    return
end
[Bi, ok] = inv3(Bm);
if ~ok
    F = intval(zeros(2));
    return
end
Z = (Cm*Bi*Cm' + (Cm*Bi*Cm')')/2;
F = Xfull - nu*Y - t*Z;
F = (F + F')/2;
end
