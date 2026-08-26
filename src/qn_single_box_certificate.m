function results = qn_single_box_certificate(face_box, radial_grid, t_levels)
% Certified INTLAB second-order difference-quotient test on one face box.
% Implements Algorithm 1, "Certificate on B," in Appendix B.4.
%
% Current paper formulas:
%   λᵢ(te) = π² + tνᵢ(t,e), and L is the scaled cluster sum.           (22)
%   P(t,e)=ν₁(t,e)ν₂(t,e) has the continuous extension in (23).
%   S(t,e) := −π²L − 2P + 2(eₐ²+e_d²){π⁴+t²(π²L+P)}.                  (24)
%   C_t(ν) is the integral difference-quotient pencil.                (28)
%   Φ_t(ν,e) and F_t(ν) define the 2×2 Schur complement.              (29)
%   det F_t(ν)=0 is equivalent to the full generalized eigenproblem.  (31)
%   F̂_t(ν) and Ψ_t(ν,e) are the affine model and remainder.           (32)
%   d̃₁ is defined in (33), and det F̂_t is written in (34).
%   L and P are enclosed by (38), with the refinement in (40)-(41).
%
% Code variable map:
%   face_box                  -> chart box from χ_{r,σ} in (42)
%   t                         -> radial interval T = [t₋,t₊]
%   I_nu, I_nu_tight          -> candidate interval [ν₋,ν₊]
%   d2_B, d1_over_t_B, d0_B   -> d₂, d̃₁, d₀ in (33)-(34)
%   S_core_B, S_padding       -> central part and error radius of [S]_B
%                                obtained from (24) and (38)
%   S_B                       -> final [S]_B enclosure in Algorithm 1
%   E_B                       -> root-error radius tη/β used after (39)
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
PI2 = FR.PI2; V5 = FR.V5;
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

% Direction enclosures for B = (E intersect S^3) x T, Eq. (42).
EI = dq2_face_direction(xI, axisdim, sgn);   % e-box E_j, evaluated as INTLAB intervals
xhp = hessianinit(xc);
EP = dq2_face_direction(xhp, axisdim, sgn);  % e at center, evaluated as INTLAB Hessian variables
xhb = hessianinit(xI);
EG = dq2_face_direction(xhb, axisdim, sgn);  % e over E_j, Hessian enclosure for centered form

% Taylor-quotient remainder bounds implement the integral in C_t from (28),
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

ad2I = sqr(EI(1)) + sqr(EI(4));               % e_a^2+e_d^2 in (24), not the area q(p)
ad2P = EP(1)*EP(1) + EP(4)*EP(4);
ad2G = EG(1)*EG(1) + EG(4)*EG(4);

% X0 structural quantities (interval chain)
KM1 = CKI{2} - PI2*CMI{2};                  % D(K-pi^2 M)(0)[e], symbolic coefficients -> interval
X0I = V5'*KM1*V5;                            % E_0 block of the first-order pencil

results = repmat(struct('ok',true,'infS',inf,'reason','skip','detail','','lam1',[],'lam2',[], ...
                        'lam3',[],'mass_lower',-inf,'lambda3_gap_lower',-inf, ...
                        'root_identity_used',false,'coarse_infS',-inf, ...
                        'root_identity_gain',0,'veigs_indices',[]), m, 1);
for k = t_levels(:)'
    results(k).ok = false; results(k).infS = -inf; results(k).reason = '';
    t = dq2_interval_hull(radial_edges(k), radial_edges(k+1));
    % -- Single-box certificate: interval blocks and Schur data -----------
    [d1_over_t_I, d2_I, d0_I, S_core_I, schur] = qn_single_box_quantities(CKI, CMI, CbI, EK, EM, Eb, ad2I, t, FR); % whole-box interval evaluation
    if isempty(schur), results(k).reason = 'affine_perp_inverse'; continue; end
    Mte = M0 + t*schur.DtM;                   % M(te), positive-definite condition in (36)
    Kte = K0 + t*schur.DtK;                   % K(te), from the same Taylor quotient
    mass_lower = gershgorin_lower_bound(Mte);
    if ~(isfinite(mass_lower) && mass_lower > 0)
        results(k).reason = 'mass_not_pd'; continue;
    end
    % ---------- centered chains (hessian at center / over box) ----------
    [d1_over_t_P, d2_P, d0_P, S_core_P] = qn_single_box_quantities(CKP, CMP, CbP, EK, EM, Eb, ad2P, t, FR); % center Hessian object
    [d1_over_t_G, d2_G, d0_G, S_core_G] = qn_single_box_quantities(CKG, CMG, CbG, EK, EM, Eb, ad2G, t, FR); % box Hessian enclosure
    if isempty(d1_over_t_P) || isempty(d1_over_t_G), results(k).reason = 'centered_inverse'; continue; end
    dx = xI - xc;
    d1_over_t_B = intersect(centered_hessian_interval(d1_over_t_P, d1_over_t_G, dx), d1_over_t_I); % d̃₁, (33)-(34)
    d2_B = intersect(centered_hessian_interval(d2_P, d2_G, dx), d2_I);                             % d₂, (34)
    d0_B = intersect(centered_hessian_interval(d0_P, d0_G, dx), d0_I);                             % d₀, (34)
    S_core_B = intersect(centered_hessian_interval(S_core_P, S_core_G, dx), S_core_I);             % central [S]_B part from (24), (38)
    if isnan(inf(d1_over_t_B)) || isnan(inf(d2_B)) || isnan(inf(d0_B)) || isnan(inf(S_core_B))
        results(k).reason = 'centered_intersection'; continue;
    end
    if inf(d2_B) <= 0, results(k).reason = 'affine_root_leading_coefficient'; continue; end
    % Roots of det Fhat_t(mu) = d2*mu^2 - t*dtilde1*mu + d0, Eq. (34).
    d1_B = t*d1_over_t_B;
    DiscC = max(sqr(d1_B) - 4*d2_B*d0_B, intval(0));
    rt = sqrt(DiscC);
    mu1 = (d1_B - rt)/(2*d2_B); mu2 = (d1_B + rt)/(2*d2_B);
    beta = inf(lammin2(schur.Bbar));          % beta in (37), interval lower bound
    if beta <= 0, results(k).reason = 'affine_slope_not_pd'; continue; end
    pad = upper_initial_pad(schur.al, schur.be); % adaptive window margin
    ok = true; reason = ''; eta = inf; E_B = inf; I_nu = [];
    nw = 6;   % window subdivisions for (C2)/(C5)/kappa2
    L1 = -schur.Nm + schur.C0*schur.B0i*schur.MW;
    nrmL1 = spnorm(L1);
    window_closed = false;
    for it = 1:3
        I_nu = padded_hull(mu1, mu2, pad);
        [kappa2, ok, reason] = schur_residual_norm_bound(schur, t, I_nu, nw, nrmL1); % interval bound for ||F_t - Fhat_t||/t
        if ~ok, break; end
        eta = upper_eta(I_nu, t, kappa2);       % eta in (37)
        E_B = upper_ebar(t, eta, beta);       % root-error radius t*eta/beta after (39)
        if upper_plus_margin(E_B) <= pad
            window_closed = true;
            break
        end
        pad = upper_scaled_margin(E_B);
    end
    if ~ok, results(k).reason = reason; continue; end
    if ~window_closed, results(k).reason = 'root_window'; continue; end
    % bootstrap: roots now known to lie in hull(mu1,mu2) +/- E_B;
    % re-evaluate the residual constants on that tighter range
    for bs = 1:2
        I_nu_tight = intersect(padded_hull(mu1, mu2, E_B), I_nu);
        [kapr, ok, reason] = schur_residual_norm_bound(schur, t, I_nu_tight, nw, nrmL1);
        if ~ok, break; end
        eta_n = upper_eta(I_nu_tight, t, kapr);
        ebar_n = upper_ebar(t, eta_n, beta);
        if ebar_n >= E_B, break; end
        eta = eta_n; E_B = ebar_n;
    end
    if ~ok, results(k).reason = reason; continue; end
    ok = true;
    Xfull = schur.X0 + t*schur.Xt;            % E_0 block of the exact Schur pencil, interval
    wedges = interval_edges(I_nu, nw);
    for w = 1:nw
        Wj = dq2_interval_hull(wedges(w), wedges(w+1));
        BWj = schur.B0 - (t*Wj)*schur.MW;     % D_perp+t*C_{t,perp,perp}(nu), (36)
        if ~posdef3(BWj), ok = false; reason = 'perp_block_not_pd'; break; end
        [BWji, okinv] = inv3(BWj);
        if ~okinv, ok = false; reason = 'perp_block_inverse'; break; end
        CWj = schur.C0 - (t*Wj)*schur.Nm;
        dZ = -(schur.Nm*BWji*CWj' + CWj*BWji*schur.Nm') + CWj*BWji*schur.MW*BWji*CWj';
        mono = schur.Y + t*symm(dZ);          % -partial_nu F_t(nu), (36)
        if inf(lammin2(symm(mono))) <= 0, ok = false; reason = 'schur_derivative'; break; end
    end
    if ~ok, results(k).reason = reason; continue; end
    F_at_nu_minus = schur_matrix_at_nu(Xfull, schur.Y, schur.C0, schur.B0, schur.Nm, schur.MW, t, intval(inf(I_nu))); % F_t(nu_-), (36)
    F_at_nu_plus = schur_matrix_at_nu(Xfull, schur.Y, schur.C0, schur.B0, schur.Nm, schur.MW, t, intval(sup(I_nu))); % F_t(nu_+), (36)
    if ~(inf(F_at_nu_minus(1,1)) > 0 && inf(det2(F_at_nu_minus)) > 0), results(k).reason = 'left_schur_not_pd'; continue; end
    if ~(sup(F_at_nu_plus(1,1)) < 0 && inf(det2(F_at_nu_plus)) > 0), results(k).reason = 'right_schur_not_nd'; continue; end
    sbar = upper_abs(I_nu);
    L_padding = sup(intval(2)*intval(eta)/intval(beta)); % L enclosure radius in (38)
    P_padding = sup(intval(2)*intval(sbar)*intval(E_B) + sqr(intval(E_B))); % product enclosure radius in (38)
    S_padding = sup(intval(L_padding)*PI2*(intval(1) + intval(2)*sqr(t)) + ...
               intval(2)*intval(P_padding)*(intval(1) + sqr(t)));
    S_B = S_core_B + midrad(0, S_padding);    % [S]_B from (24) and (38)
    nu1_interval = mu1 + midrad(0, E_B);      % ν₁ enclosure from the certified root window
    nu2_interval = mu2 + midrad(0, E_B);      % ν₂ enclosure from the certified root window
    % Save the coarse eta/beta enclosure before applying the exact root
    % identity below.  The identity is part of the certified box algorithm;
    % keeping both bounds makes its contribution auditable in diagnostics.
    S_B_coarse = S_B;
    root_identity_used = false;
    % ---- signed second-order correction via the exact 2x2 root identity ----
    gap = nu2_interval - nu1_interval;
    if inf(gap) > 0
        Apen = schur.X0 + t*schur.G;
        root_shift_over_t = intval(zeros(2,1)); root_identity_ok = true;
        nu_intervals = {nu1_interval, nu2_interval};
        for kk = 1:2
            nu_k = nu_intervals{kk};
            B_at_nu_k = schur.B0 - (t*nu_k)*schur.MW;
            [B_at_nu_k_inv, okinv] = inv3(B_at_nu_k);
            if ~okinv, root_identity_ok = false; break; end
            L2_at_nu_k = -B_at_nu_k_inv*schur.Nm' + schur.B0i*schur.MW*B_at_nu_k_inv*schur.C0';
            L1L2_at_nu_k = L1*L2_at_nu_k;
            L1L2_at_nu_k = (L1L2_at_nu_k + L1L2_at_nu_k')/2;
            residual_over_t = sqr(t)*sqr(nu_k)*L1L2_at_nu_k;       % E(ν)/t with E = t³ν²L₁L₂
            affine_at_nu_k = Apen - nu_k*schur.Bbar;
            adj_affine = [affine_at_nu_k(2,2), -affine_at_nu_k(1,2); -affine_at_nu_k(2,1), affine_at_nu_k(1,1)];
            trace_adj_residual = adj_affine(1,1)*residual_over_t(1,1) + adj_affine(1,2)*residual_over_t(2,1) + ...
                                 adj_affine(2,1)*residual_over_t(1,2) + adj_affine(2,2)*residual_over_t(2,2);
            det_residual_over_t = t*sqr(sqr(t))*sqr(sqr(nu_k))*det2(L1L2_at_nu_k);   % det(E)/t = t⁵ν⁴det(L₁L₂)
            if kk == 1, root_gap = nu1_interval - nu2_interval; else, root_gap = nu2_interval - nu1_interval; end
            root_shift_over_t(kk) = (trace_adj_residual - det_residual_over_t)/(d2_B*root_gap);
        end
        if root_identity_ok
            root_identity_used = true;
            L_correction = root_shift_over_t(1) + root_shift_over_t(2); % (δ₁+δ₂)/t
            P_correction = t*(nu1_interval*root_shift_over_t(2) + nu2_interval*root_shift_over_t(1)) + ...
                           sqr(t)*root_shift_over_t(1)*root_shift_over_t(2);
            S_correction = -PI2*L_correction - 2*P_correction + 2*ad2I*sqr(t)*(PI2*L_correction + P_correction); % (41) substituted in S from (24)
            S_B_from_root_identity = S_core_B + S_correction;
            S_B = intersect(S_B, S_B_from_root_identity);
            nu1_interval = intersect(nu1_interval, mu1 + t*root_shift_over_t(1));
            nu2_interval = intersect(nu2_interval, mu2 + t*root_shift_over_t(2));
        end
    end
    if isnan(inf(nu1_interval)) || isnan(sup(nu1_interval)) || ...
            isnan(inf(nu2_interval)) || isnan(sup(nu2_interval))
        results(k).reason = 'root_intersection'; continue;
    end
    lam1 = PI2 + t*nu1_interval;              % refined λ₁(te), (22)
    lam2 = PI2 + t*nu2_interval;              % refined λ₂(te), (22)
    % Independently target indices 1 and 3 of the complete 5x5 interval
    % pencil with verified veigs/veig calls. The DQ2 Schur analysis remains necessary
    % for the double cluster and S(t,e); acceptance additionally requires
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
    results(k).root_identity_used = root_identity_used;
    results(k).coarse_infS = inf(S_B_coarse);
    results(k).root_identity_gain = inf(S_B) - inf(S_B_coarse);
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

function W = padded_hull(mu1, mu2, pad)
% [inf(mu1)-pad, sup(mu2)+pad], with endpoint arithmetic done by INTLAB.
lo = inf(intval(inf(mu1)) - intval(pad));
hi = sup(intval(sup(mu2)) + intval(pad));
W = infsup(lo, hi);
end

function eta = upper_eta(W, t, kappa2)
% Upward-rounded eta bound for (37).
sbar = upper_abs(W);
eta = sup(sqr(intval(sbar)) * sqr(intval(t)) * intval(kappa2));
end

function ebar = upper_ebar(t, eta, beta)
% Upward-rounded root-error radius t*eta/beta used after (39).
ebar = sup(intval(t) * intval(eta) / intval(beta));
end

function y = upper_plus_margin(x)
y = sup(intval(x) + intval('1e-6'));
end

function y = upper_scaled_margin(x)
y = sup(intval('1.5') * (intval(x) + intval('1e-6')));
end

function n = fronorm(A)
n = sup(sqrt(sum(sum(sqr(intval(A))))));
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
% kappa2 = ||L1|| * sup over I_nu of ||L2(mu)||, by window subdivision
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

function P = schur_matrix_at_nu(Xfull, Y, C0, B0, Nm, MW, t, mu)
Cm = C0 - (t*mu)*Nm;
Bm = B0 - (t*mu)*MW;
[Bi, ok] = inv3(Bm);
if ~ok
    P = intval([-1 0; 0 -1]);
    return
end
Z = (Cm*Bi*Cm' + (Cm*Bi*Cm')')/2;
P = Xfull - mu*Y - t*Z;
P = (P + P')/2;
end
