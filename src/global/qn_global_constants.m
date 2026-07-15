function C = qn_global_constants()
% Paper-aligned constants for the certified global cover.
% Every non-integer proof constant is parsed by INTLAB from a decimal string.

PI = intval('pi');
C.pi2 = PI^2;
% Explicit decimal upper bound for the first positive zero of J_0.
C.j01_upper = intval('2.404825557695774');
C.T = intval('4') * C.j01_upper^2 / C.pi2;
C.lambda_min_H2 = intval('3232') / (intval('27') * PI^6);
C.rho_local = C.lambda_min_H2;
C.rho_seam = C.lambda_min_H2 / intval('2');
C.ulp_pad = intval('1e-12');
C.gap_threshold = 1.5; % non-certified frame-selection heuristic
C.max_depth = 60;

s = sqrt(C.T);
C.pc_half = [sqrt(C.T / intval('2') - intval('1')); ...
    (s - intval('2') + sqrt(s^2 + intval('4')*s - intval('4'))) / intval('2'); ...
    (s - intval('2') + sqrt(s^2 + intval('4')*s - intval('4'))) / intval('2'); ...
    sqrt(C.T / intval('2')) - intval('1')];
end
