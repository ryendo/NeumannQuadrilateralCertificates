function results = dq2_retest_sign_box(child_box, radial_grid, t_levels, parent_data)
% Cheap re-test of S-only failures on a sub-box.
% This only rechecks the final sign condition inf S_B > 0 from
% Algorithm 1, Steps 13-15; the ancestor box already certified the
% Schur/root conditions.
%
% PDF notation:
%   S(t,e) and S_B are Eqs. (13), (67).  Their coefficients come from the
%   same Vieta data d_2, tilde_d_1, d_0 in Eqs. (59)-(60).
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

PI = intval('pi');
PI2 = sqr(PI);
TC = dq2_load_taylor_coefficients();
SQ2 = sqrt(intval(2));
V5 = intval(zeros(5,2)); V5(1,1) = SQ2; V5(2,2) = SQ2;
W5 = intval(zeros(5,3)); W5(3,1) = intval(2); W5(4,2) = SQ2; W5(5,3) = SQ2;
D0 = intval(zeros(3)); D0(1,1) = PI2; D0(2,2) = 3*PI2; D0(3,3) = 3*PI2;
I2 = intval(eye(2)); I3 = intval(eye(3));
FR = struct('V5',V5,'W5',W5,'D0',D0,'I2',I2,'I3',I3,'PI2',PI2);

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
EP = e_from_x(xhp, axisdim, sgn);
EK = midrad(zeros(5), parent_data.REM.BK); EM = midrad(zeros(5), parent_data.REM.BM);
Eb = midrad(zeros(5,1), parent_data.REM.Bb);
[CKP, CMP, CbP] = dq2_evaluate_taylor_coefficients_vectorized(TC, EP(1), EP(2), EP(3), EP(4)); % symbolic coefficients -> center Hessian data
qP = EP(1)*EP(1) + EP(4)*EP(4);

m = size(radial_grid, 2) - 1;
results = repmat(struct('ok',true,'infS',inf,'reason','skip'), m, 1);
for k = t_levels(:)'
    t = interval_hull(radial_grid(k), radial_grid(k+1));
    [~, d2_P, ~, S_core_center_data] = algorithm1_scalars_sign_only(CKP, CMP, CbP, EK, EM, Eb, qP, t, FR); % centered S_B data, Eq. (67)
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
    S_B_child = S_core_child + midrad(0, parent_data.S_padding(k)); % final sign test, Eq. (67)
    results(k).infS = inf(S_B_child);
    results(k).ok = inf(S_B_child) > 0;
    if ~results(k).ok, results(k).reason = 'S'; end
end
end

function X = interval_hull(a, b)
% Hull of two scalar endpoints, accepting either double or intval inputs.
ia = intval(a); ib = intval(b);
X = infsup(inf(ia), sup(ib));
end

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

function [d1_over_t, d2, d0, S_core] = algorithm1_scalars_sign_only(CK, CM, Cb, EK, EM, Eb, q, t, FR)
% Same as dq2_algorithm1_box>algorithm1_scalars, without block output.
% The coefficients are those of Eq. (59), with d1 = t*d1_over_t from Eq. (60).
PI2 = FR.PI2; V5 = FR.V5; W5 = FR.W5; D0 = FR.D0; I2 = FR.I2; I3 = FR.I3;
J = 9;
RK = CK{3}; RMv = CM{3}; Rb = Cb{3};
tp = t;
for j = 3:J
    RK = RK + tp*CK{j+1}; RMv = RMv + tp*CM{j+1}; Rb = Rb + tp*Cb{j+1};
    tp = tp*t;
end
RK = RK + tp*EK; RMv = RMv + tp*EM; Rb = Rb + tp*Eb;
den = 1 - sqr(t)*q;                            % area factor |Q_te|, Eq. (9)
wv = Cb{2} + t*Rb;
RM = RMv - (wv*wv')/den;
Xt = V5'*(RK - PI2*RM)*V5;
DtM = CM{2} + t*RM; DtK = CK{2} + t*RK;
Th = DtK - PI2*DtM;
Yd = V5'*DtM*V5;
Nm = V5'*DtM*W5;
MW = I3 + t*(W5'*DtM*W5);
C0 = V5'*Th*W5;
B0 = D0 + t*(W5'*Th*W5);
det_B0 = B0(1,1)*(B0(2,2)*B0(3,3)-B0(2,3)*B0(3,2)) ...
  - B0(1,2)*(B0(2,1)*B0(3,3)-B0(2,3)*B0(3,1)) ...
  + B0(1,3)*(B0(2,1)*B0(3,2)-B0(2,2)*B0(3,1));
det_B0_interval = det_B0.x;
if ~((inf(det_B0_interval) > 0) || (sup(det_B0_interval) < 0))
    d1_over_t = []; d2 = []; d0 = []; S_core = [];
    return
end
B0i = [ B0(2,2)*B0(3,3)-B0(2,3)*B0(3,2), B0(1,3)*B0(3,2)-B0(1,2)*B0(3,3), B0(1,2)*B0(2,3)-B0(1,3)*B0(2,2);
        B0(2,3)*B0(3,1)-B0(2,1)*B0(3,3), B0(1,1)*B0(3,3)-B0(1,3)*B0(3,1), B0(1,3)*B0(2,1)-B0(1,1)*B0(2,3);
        B0(2,1)*B0(3,2)-B0(2,2)*B0(3,1), B0(1,2)*B0(3,1)-B0(1,1)*B0(3,2), B0(1,1)*B0(2,2)-B0(1,2)*B0(2,1) ]/det_B0;
Z0 = (C0*B0i*C0' + (C0*B0i*C0')')/2;
CB0 = C0*B0i;
Zb1_unsym = -(Nm*B0i*C0' + C0*B0i*Nm') + CB0*MW*B0i*C0';
Zb1 = (Zb1_unsym + Zb1_unsym')/2;
X0 = V5'*(CK{2} - PI2*CM{2})*V5;
al = (X0(1,1) - X0(2,2))/2;
be = (X0(1,2) + X0(2,1))/2;
G = Xt - Z0;
CBm = Yd + t*Zb1;
BB = I2 + t*CBm;
d2 = BB(1,1)*BB(2,2) - BB(1,2)*BB(2,1);
d1_over_t = (G(1,1)+G(2,2)) + al*(CBm(2,2)-CBm(1,1)) - 2*be*CBm(1,2) ...
      + t*(G(1,1)*CBm(2,2) + G(2,2)*CBm(1,1) - 2*G(1,2)*CBm(1,2));
d0 = -(al*al+be*be) + t*(al*(G(2,2)-G(1,1)) - 2*be*G(1,2)) ...
     + sqr(t)*(G(1,1)*G(2,2) - G(1,2)*G(1,2));
PI4 = sqr(PI2);
S_core = -(PI2*d1_over_t + 2*d0)/d2 + 2*q*(PI4 + sqr(t)*(PI2*d1_over_t + d0)/d2); % S_B core, Eq. (67)
end
