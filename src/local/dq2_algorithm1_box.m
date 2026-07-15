function [results, REM] = dq2_algorithm1_box(face_box, m, radial_grid, parent_remainder_bounds, t_levels)
% Certified INTLAB second-order difference-quotient test on one face box.
% This implements Algorithm 1 (Single-box certificate) from Sec. 4.1.
%
% PDF formulas from DQ2-quadrilateral-note-ver4.pdf:
%   K(p)x = λM(p)x.                                                   (3)
%   λᵢ(te) = π² + tνᵢ(t,e).                                           (5)
%   L(t,e) = (ν₁(t,e)+ν₂(t,e))/t,   P(t,e)=ν₁(t,e)ν₂(t,e).             (6)
%   q(e)=eₐ²+e_d²,   Π_t=λ₁(te)λ₂(te).                                (8)
%   S(t,e) := −π²L(t,e) − 2P(t,e) + 2q(e)Π_t.                         (13)
%   Δ_tK(e), Δ_tM(e) use the limiting definitions at t=0.        (29)-(30)
%   C_t(ν) := Δ_tK(e) − π²Δ_tM(e) − νM(te).                           (32)
%   F_t(ν) is the 2×2 Schur complement.                               (39)
%   F̂_t(ν) = F_t(0) + ν∂_νF_t(0).                                    (44),(58)
%   det F̂_t(ν) = d₂ν² − d₁ν + d₀,   d₁ = t d̃₁.                (45),(59),(60)
%   L_B, P_B, q_B, Π_B, S_B are Algorithm 1 interval enclosures.  (63)-(67)
%
% Code variable map:
%   face_box                  -> chart box E in B = (E ∩ S³) × T, Eq. (42)
%   t                         -> radial interval T = [t₋,t₊]
%   I_nu, I_nu_tight                 -> candidate ν intervals Iν = [ν₋,ν₊]
%   d2_B, d1_over_t_B, d0_B   -> d₂, d̃₁, d₀ in Eqs. (59)-(60)
%   S_core_B                  -> unpadded core of S_B in Eq. (67)
%   S_padding                 -> interval radius added to S_core_B in Eq. (67)
%   L_B, P_B, S_B             -> interval enclosures in Eqs. (63),(64),(67)
%   E_B                       -> E_B = Tη/β in Eq. (62)
%
% Evaluation:
%   All rigorous quantities are INTLAB intervals.  Taylor coefficients of
%   K, M, and b are precomputed symbolically and loaded by
%   dq2_load_taylor_coefficients.m; they are then evaluated either as INTLAB
%   intervals (whole box) or as INTLAB Hessian objects (centered form).
%
% face_box = [axisdim sgn lo1 hi1 lo2 hi2 lo3 hi3].
% parent_remainder_bounds: optional bounds inherited from a parent box
% (valid for any sub-box). t_levels: optional subset of radial intervals.

PI = intval('pi');
PI2 = sqr(PI); PI4 = sqr(PI2);
rho = 3232/(27*PI2*PI2*PI2);                  % rho# in Sec. 1
tmax = intval(sup(rho))*intval('1.000000000001');
TC = dq2_load_taylor_coefficients();
SQ2 = sqrt(intval(2));
V5 = intval(zeros(5,2)); V5(1,1) = SQ2; V5(2,2) = SQ2;          % E_0 basis, Eq. (18), interval constants
W5 = intval(zeros(5,3)); W5(3,1) = intval(2); W5(4,2) = SQ2; W5(5,3) = SQ2; % E_perp basis, Eq. (18)
D0 = intval(zeros(3)); D0(1,1) = PI2; D0(2,2) = 3*PI2; D0(3,3) = 3*PI2;    % D_perp, Eq. (22)
I2 = intval(eye(2)); I3 = intval(eye(3));
M0 = intval(diag([1/2,1/2,1/4,1/2,1/2]));
FR = struct('V5',V5,'W5',W5,'D0',D0,'I2',I2,'I3',I3,'PI2',PI2);

axisdim = face_box(1); sgn = face_box(2);
xI = intval(zeros(3,1));
for d = 1:3
    xI(d) = infsup(face_box(1+2*d), face_box(2+2*d));
end
xc = intval(mid(xI));

if nargin < 3 || isempty(radial_grid)
    radial_edges = tmax * (intval(0:m)/intval(m));
else
    radial_edges = radial_grid; m = size(radial_edges, 2) - 1;
end

% Direction enclosures for B = (E intersect S^3) x T, Eq. (42).
EI = e_from_x(xI, axisdim, sgn);            % e-box E_j, evaluated as INTLAB intervals
xhp = hessianinit(xc);
EP = e_from_x(xhp, axisdim, sgn);           % e at center, evaluated as INTLAB Hessian variables
xhb = hessianinit(xI);
EG = e_from_x(xhb, axisdim, sgn);           % e over E_j, Hessian enclosure for centered form

% Taylor-quotient remainder bounds implement Delta_t K and Delta_t M in
% Eqs. (29)-(30), avoiding division by an interval containing t = 0.
if nargin >= 4 && ~isempty(parent_remainder_bounds)
    BKr = parent_remainder_bounds.BK; BMr = parent_remainder_bounds.BM; Bbr = parent_remainder_bounds.Bb; % inherited interval remainders R_K, R_M, R_b
else
    [BKr, BMr, Bbr] = dq2_bound_taylor_remainder(EI, radial_edges(end), 2, 10); % direct INTLAB Taylor remainder
end
REM = struct('BK', BKr, 'BM', BMr, 'Bb', Bbr);
EK = midrad(zeros(5), BKr); EM = midrad(zeros(5), BMr); Eb = midrad(zeros(5,1), Bbr); % interval R_K, R_M, R_b
if nargin < 5 || isempty(t_levels)
    t_levels = 1:m;
end

% Symbolic Taylor coefficient table evaluated in three rigorous modes:
[CKI, CMI, CbI] = dq2_evaluate_taylor_coefficients(TC, EI(1), EI(2), EI(3), EI(4)); % whole-box INTLAB intervals
[CKP, CMP, CbP] = dq2_evaluate_taylor_coefficients(TC, EP(1), EP(2), EP(3), EP(4)); % center value/gradient
[CKG, CMG, CbG] = dq2_evaluate_taylor_coefficients(TC, EG(1), EG(2), EG(3), EG(4)); % Hessian enclosure over box

qI = sqr(EI(1)) + sqr(EI(4));                 % q(e) in Eq. (8)
qP = EP(1)*EP(1) + EP(4)*EP(4);
qG = EG(1)*EG(1) + EG(4)*EG(4);

% X0 structural quantities (interval chain)
KM1 = CKI{2} - PI2*CMI{2};                  % D(K-pi^2 M)(0)[e], symbolic coefficients -> interval
X0I = V5'*KM1*V5;                            % E_0 block of the first-order pencil

results = repmat(struct('ok',true,'infS',inf,'reason','skip','lam1',[],'lam2',[], ...
                        'S_interval',[],'L_interval',[],'E_B',0,'nu_window',[], ...
                        'S_padding',0,'S_core_hessian',[], ...
                        'root_identity_used',false,'algorithm1_infS',-inf, ...
                        'root_identity_gain',0), m, 1);
for k = t_levels(:)'
    results(k).ok = false; results(k).infS = -inf; results(k).reason = '';
    t = interval_hull(radial_edges(k), radial_edges(k+1));
    % ---------- Algorithm 1 Steps 1-3: interval blocks and Schur data ----
    [d1_over_t_I, d2_I, d0_I, S_core_I, schur] = algorithm1_scalars(CKI, CMI, CbI, EK, EM, Eb, qI, t, FR); % whole-box interval evaluation
    if isempty(schur), results(k).reason = 'B0inv'; continue; end
    % ---------- centered chains (hessian at center / over box) ----------
    [d1_over_t_P, d2_P, d0_P, S_core_P] = algorithm1_scalars(CKP, CMP, CbP, EK, EM, Eb, qP, t, FR); % center Hessian object
    [d1_over_t_G, d2_G, d0_G, S_core_G] = algorithm1_scalars(CKG, CMG, CbG, EK, EM, Eb, qG, t, FR); % box Hessian enclosure
    if isempty(d1_over_t_P) || isempty(d1_over_t_G), results(k).reason = 'cinv'; continue; end
    dx = xI - xc;
    d1_over_t_B = intersect(centered_hessian_interval(d1_over_t_P, d1_over_t_G, dx), d1_over_t_I); % d̃₁, Eq. (60)
    d2_B = intersect(centered_hessian_interval(d2_P, d2_G, dx), d2_I);                             % d₂, Eq. (59)
    d0_B = intersect(centered_hessian_interval(d0_P, d0_G, dx), d0_I);                             % d₀, Eq. (59)
    S_core_B = intersect(centered_hessian_interval(S_core_P, S_core_G, dx), S_core_I);             % unpadded S_B, Eq. (67)
    if isnan(inf(d1_over_t_B)) || isnan(inf(d2_B)) || isnan(inf(d0_B)) || isnan(inf(S_core_B))
        results(k).reason = 'isect'; continue;
    end
    if inf(d2_B) <= 0, results(k).reason = 'c2'; continue; end
    % Roots of det Fhat_t(mu) = d2*mu^2 - d1*mu + d0, Eqs. (59)-(60).
    d1_B = t*d1_over_t_B;
    DiscC = max(sqr(d1_B) - 4*d2_B*d0_B, intval(0));
    rt = sqrt(DiscC);
    mu1 = (d1_B - rt)/(2*d2_B); mu2 = (d1_B + rt)/(2*d2_B);
    beta = inf(lammin2(schur.Bbar));          % beta in Eq. (61), interval lower bound
    if beta <= 0, results(k).reason = 'beta'; continue; end
    pad = upper_initial_pad(schur.al, schur.be); % adaptive window margin
    ok = true; reason = ''; eta = inf; E_B = inf; I_nu = [];
    nw = 6;   % window subdivisions for (C2)/(C5)/kappa2
    L1 = -schur.Nm + schur.C0*schur.B0i*schur.MW;
    nrmL1 = spnorm(L1);
    for it = 1:3
        I_nu = padded_hull(mu1, mu2, pad);
        [kappa2, ok] = schur_residual_norm_bound(schur, t, I_nu, nw, nrmL1); % interval bound for ||F_t - Fhat_t||/t
        if ~ok, reason = 'BWinv'; break; end
        eta = upper_eta(I_nu, t, kappa2);       % eta for Eq. (61)
        E_B = upper_ebar(t, eta, beta);       % E_B in Eq. (62)
        if upper_plus_margin(E_B) <= pad, break; end
        pad = upper_scaled_margin(E_B);
    end
    if ~ok, results(k).reason = reason; continue; end
    % bootstrap: roots now known to lie in hull(mu1,mu2) +/- E_B;
    % re-evaluate the residual constants on that tighter range
    for bs = 1:2
        I_nu_tight = intersect(padded_hull(mu1, mu2, E_B), I_nu);
        [kapr, ok] = schur_residual_norm_bound(schur, t, I_nu_tight, nw, nrmL1);
        if ~ok, break; end
        eta_n = upper_eta(I_nu_tight, t, kapr);
        ebar_n = upper_ebar(t, eta_n, beta);
        if ebar_n >= E_B, break; end
        eta = eta_n; E_B = ebar_n;
    end
    ok = true;
    Xfull = schur.X0 + t*schur.Xt;            % E_0 block of the exact Schur pencil, interval
    wedges = interval_edges(I_nu, nw);
    for w = 1:nw
        Wj = interval_hull(wedges(w), wedges(w+1));
        BWj = schur.B0 - (t*Wj)*schur.MW;     % A_t^{perp,perp}(nu), Eq. (57)
        if ~posdef3(BWj), ok = false; reason = 'W2'; break; end
        [BWji, okinv] = inv3(BWj);
        if ~okinv, ok = false; reason = 'BWinv'; break; end
        CWj = schur.C0 - (t*Wj)*schur.Nm;
        dZ = -(schur.Nm*BWji*CWj' + CWj*BWji*schur.Nm') + CWj*BWji*schur.MW*BWji*CWj';
        mono = schur.Y + t*symm(dZ);          % -partial_nu F_t(nu), Eq. (57)
        if inf(lammin2(symm(mono))) <= 0, ok = false; reason = 'W5'; break; end
    end
    if ~ok, results(k).reason = reason; continue; end
    F_at_nu_minus = schur_matrix_at_nu(Xfull, schur.Y, schur.C0, schur.B0, schur.Nm, schur.MW, t, intval(inf(I_nu))); % F_t(nu_-), Eq. (57)
    F_at_nu_plus = schur_matrix_at_nu(Xfull, schur.Y, schur.C0, schur.B0, schur.Nm, schur.MW, t, intval(sup(I_nu))); % F_t(nu_+), Eq. (57)
    if ~(inf(F_at_nu_minus(1,1)) > 0 && inf(det2(F_at_nu_minus)) > 0), results(k).reason = 'W3'; continue; end
    if ~(sup(F_at_nu_plus(1,1)) < 0 && inf(det2(F_at_nu_plus)) > 0), results(k).reason = 'W4'; continue; end
    Mte = M0 + t*schur.DtM;                   % M(te), Eqs. (3), (57), positive-definite check
    gersh = true;
    for i = 1:5
        jj = [1:i-1, i+1:5];
        % mag(...) gives componentwise floating-point upper bounds.  Convert
        % them back to intervals before summing so the row-radius sum is
        % guaranteed to be rounded upward.
        off = sum(intval(mag(Mte(i,jj))));
        if inf(Mte(i,i) - off) <= 0, gersh = false; break; end
    end
    if ~gersh, results(k).reason = 'W1'; continue; end
    sbar = upper_abs(I_nu);
    L_padding = sup(intval(2)*intval(eta)/intval(beta)); % L_B padding, Eq. (63)
    P_padding = sup(intval(2)*intval(sbar)*intval(E_B) + sqr(intval(E_B))); % P_B padding, Eq. (64)
    S_padding = sup(intval(L_padding)*PI2*(intval(1) + intval(2)*sqr(t)) + ...
               intval(2)*intval(P_padding)*(intval(1) + sqr(t)));
    S_B = S_core_B + midrad(0, S_padding);    % S_B interval, Eq. (67)
    L_B = d1_over_t_B/d2_B + midrad(0, L_padding); % L_B interval, Eq. (63)
    P_B = d0_B/d2_B + midrad(0, P_padding);   % P_B interval, Eq. (64)
    nu1_interval = mu1 + midrad(0, E_B);      % ν₁ root enclosure, Eq. (53)
    nu2_interval = mu2 + midrad(0, E_B);      % ν₂ root enclosure, Eq. (53)
    % Save the coarse eta/beta enclosure before applying the exact root
    % identity below.  The identity is part of the certified box algorithm;
    % keeping both bounds makes its contribution auditable in diagnostics.
    S_B_algorithm1 = S_B;
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
            L_B_from_root_identity = d1_over_t_B/d2_B + L_correction;
            P_B_from_root_identity = d0_B/d2_B + P_correction;
            L_B = intersect(L_B, L_B_from_root_identity);
            P_B = intersect(P_B, P_B_from_root_identity);
            S_correction = -PI2*L_correction - 2*P_correction + 2*qI*sqr(t)*(PI2*L_correction + P_correction); % correction for S_B, Eq. (67)
            S_B_from_root_identity = S_core_B + S_correction;
            S_B = intersect(S_B, S_B_from_root_identity);
            S_padding = min(S_padding, mag(S_correction));
            nu1_interval = intersect(nu1_interval, mu1 + t*root_shift_over_t(1));
            nu2_interval = intersect(nu2_interval, mu2 + t*root_shift_over_t(2));
        end
    end
    lam1 = PI2 + t*nu1_interval;              % refined λ₁(te), Eq. (5)
    lam2 = PI2 + t*nu2_interval;              % refined λ₂(te), Eq. (5)
    results(k).infS = inf(S_B);
    results(k).lam1 = [inf(lam1), sup(lam1)];
    results(k).lam2 = [inf(lam2), sup(lam2)];
    results(k).S_interval = [inf(S_B), sup(S_B)];
    results(k).L_interval = [inf(L_B), sup(L_B)];
    results(k).E_B = E_B;
    results(k).nu_window = [inf(mu1), sup(mu2)];
    results(k).S_padding = S_padding;
    results(k).S_core_hessian = S_core_G;
    results(k).root_identity_used = root_identity_used;
    results(k).algorithm1_infS = inf(S_B_algorithm1);
    results(k).root_identity_gain = inf(S_B) - inf(S_B_algorithm1);
    if inf(S_B) > 0 && inf(lam1) > 0
        results(k).ok = true;
    else
        results(k).reason = 'S';
    end
end
end

% ---------- helpers ----------
function E = e_from_x(x, axisdim, sgn)
oth = setdiff(1:4, axisdim);
nrm = sqrt(1 + x(1)^2 + x(2)^2 + x(3)^2);
E = cell(4,1);
E{axisdim} = sgn/nrm;
for d = 1:3
    E{oth(d)} = x(d)/nrm;
end
E = vertcat(E{:});
end

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

function [d1_over_t, d2, d0, S_core, schur] = algorithm1_scalars(CK, CM, Cb, EK, EM, Eb, q, t, FR)
% Vieta scalars and the assembled unpadded S_B for the comparison pencil;
% class-generic (intval / hessian). schur (block data) for intval inputs.
% d2, d1_over_t, d0 are d₂, d̃₁, d₀ in Algorithm 1, Steps 6 and 9-13:
% det Fhat_t(nu) = d_2 nu^2 - d_1 nu + d_0, d_1 = t*tilde_d_1,
% Eqs. (59)-(60). S_core is Eq. (67) before adding the root-error padding.
PI2 = FR.PI2; V5 = FR.V5; W5 = FR.W5; D0 = FR.D0; I2 = FR.I2; I3 = FR.I3;
J = 9;
RK = CK{3}; RMv = CM{3}; Rb = Cb{3};          % K_2(e), M_2(e), b_2(e), symbolic coefficients
tp = t;
for j = 3:J
    RK = RK + tp*CK{j+1}; RMv = RMv + tp*CM{j+1}; Rb = Rb + tp*Cb{j+1}; % Taylor quotient for Eqs. (29)-(30)
    tp = tp*t;
end
% tp = t^8
RK = RK + tp*EK; RMv = RMv + tp*EM; Rb = Rb + tp*Eb; % interval Taylor remainders R_K, R_M, R_b
den = 1 - sqr(t)*q;                            % area factor |Q_te|, Eq. (9)
wv = Cb{2} + t*Rb;
RM = RMv - (wv*wv')/den;                      % Delta_t M(e), including rank-one correction
Xt = V5'*(RK - PI2*RM)*V5;                    % E_0 block of second-order quotient
DtM = CM{2} + t*RM; DtK = CK{2} + t*RK;       % Delta_t M(e), Delta_t K(e), Eqs. (29)-(30)
Th = DtK - PI2*DtM;                           % Delta_t K(e) - pi^2 Delta_t M(e)
Yd = V5'*DtM*V5;                              % E_0 block of Delta_t M(e)
Nm = V5'*DtM*W5;                              % E_0/E_perp block of Delta_t M(e)
MW = I3 + t*(W5'*DtM*W5);
C0 = V5'*Th*W5;                               % C_t^{0 perp}(0), Eq. (32)
B0 = D0 + t*(W5'*Th*W5);                       % A_t^{perp,perp}(0), Eqs. (34), (57)
[B0i, okinv] = inv3(B0);
if ~okinv
    d1_over_t = []; d2 = []; d0 = []; S_core = []; schur = [];
    return
end
Z0 = symm(C0*B0i*C0');                        % Schur correction in F_t, Eq. (39), interval 3x3 solve
CB0 = C0*B0i;
Zb1 = symm(-(Nm*B0i*C0' + C0*B0i*Nm') + CB0*MW*B0i*C0');
X0 = V5'*(CK{2} - PI2*CM{2})*V5;              % first-order E_0 block at t=0
al = (X0(1,1) - X0(2,2))/2;
be = (X0(1,2) + X0(2,1))/2;
G = Xt - Z0;                                  % constant part of Fhat_t
CBm = Yd + t*Zb1;                             % coefficient of -nu in Fhat_t, Eq. (58)
d2 = det2(I2 + t*CBm);                         % d₂ in Eq. (59)
d1_over_t = (G(1,1)+G(2,2)) + al*(CBm(2,2)-CBm(1,1)) - 2*be*CBm(1,2) ...
      + t*(G(1,1)*CBm(2,2) + G(2,2)*CBm(1,1) - 2*G(1,2)*CBm(1,2));
d0 = -(sqr2g(al)+sqr2g(be)) + t*(al*(G(2,2)-G(1,1)) - 2*be*G(1,2)) ...
     + sqr(t)*(G(1,1)*G(2,2) - sqr2g(G(1,2)));
PI4 = sqr(FR.PI2);
S_core = -(PI2*d1_over_t + 2*d0)/d2 + 2*q*(PI4 + sqr(t)*(PI2*d1_over_t + d0)/d2); % S_B core, Eq. (67)
if nargout > 4
    schur = struct('X0',X0,'al',al,'be',be,'Xt',Xt,'DtM',DtM,'Y',I2+t*Yd, ...
                   'Nm',Nm,'MW',MW,'C0',C0,'B0',B0,'B0i',B0i,'Bbar',I2+t*CBm, ...
                   'Z0',Z0,'G',G);
end
end

function y = sqr2g(x)
% square that works for both intval (tight sqr) and gradient
if isintval(x)
    y = sqr(x);
else
    y = x*x;
end
end

function S = symm(A)
S = (A + A')/2;
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

function X = interval_hull(a, b)
% Hull of two scalar endpoints, accepting either double or intval inputs.
ia = intval(a); ib = intval(b);
X = infsup(inf(ia), sup(ib));
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
% Upward-rounded eta bound for Eq. (61).
sbar = upper_abs(W);
eta = sup(sqr(intval(sbar)) * sqr(intval(t)) * intval(kappa2));
end

function ebar = upper_ebar(t, eta, beta)
% Upward-rounded E_B = T eta / beta from Eq. (62).
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

function [kappa2, ok] = schur_residual_norm_bound(schur, t, I_nu, nw, nrmL1)
% kappa2 = ||L1|| * sup over I_nu of ||L2(mu)||, by window subdivision
wedges = interval_edges(I_nu, nw);
kappa2 = 0; ok = true;
for w = 1:nw
    Wj = interval_hull(wedges(w), wedges(w+1));
    BWj = schur.B0 - (t*Wj)*schur.MW;
    [BWji, okinv] = inv3(BWj);
    if ~okinv, ok = false; return; end
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
