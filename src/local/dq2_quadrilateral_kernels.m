function out = dq2_quadrilateral_kernels(u, v, aa, bb, cc, dd, PI)
% Integrands of K(p), the volume part of M(p), and the rank-one mass vector
% at scalar (u,v), parameters aa..dd.
% K(p) and M(p) are the 5-by-5 matrices in the generalized eigenproblem
% K(p) x = lambda M(p) x, Eq. (3).
% Works for intval or taylor parameter types (only +,-,*,/,sin,cos used).
%
% PDF notation:
%   K(p), M(p) are Eq. (3).  The five modes are listed in Eq. (14).
%
% Evaluation:
%   This file evaluates the exact symbolic integrands directly.  Depending
%   on the input type, the same formulas produce INTLAB intervals or INTLAB
%   Taylor coefficients used to enclose Delta_t K(e), Delta_t M(e) in
%   Algorithm 1, Step 1, Eq. (54).
modes = [1 0; 0 1; 1 1; 2 0; 0 2];
X = u - aa*u - dd*v - 2*bb*(u*v);            % quadrilateral pullback x(u,v;p)
Y = v - dd*u + aa*v + 2*cc*(u*v);            % quadrilateral pullback y(u,v;p)
Xu = 1 - aa - 2*bb*v;
Xv = -(dd + 2*bb*u);
Yu = -(dd - 2*cc*v);
Yv = 1 + aa + 2*cc*u;
J = Xu*Yv - Xv*Yu;                            % Jacobian, contributes to M(p)
Lam = cell(5,1); dLu = cell(5,1); dLv = cell(5,1);
for i = 1:5
    m = modes(i,1); n = modes(i,2);
    cx = cos(m*PI*(X + 0.5)); sx = sin(m*PI*(X + 0.5));
    cy = cos(n*PI*(Y + 0.5)); sy = sin(n*PI*(Y + 0.5));
    Lam{i} = cx*cy;                           % basis function phi_i
    Lx = -m*PI*(sx*cy);
    Ly = -n*PI*(cx*sy);
    dLu{i} = Lx*Xu + Ly*Yu;
    dLv{i} = Lx*Xv + Ly*Yv;
end
out = cell(35,1);
cnt = 0;
for i = 1:5
    for j = i:5
        cnt = cnt + 1;
        num = (Xv*Xv + Yv*Yv)*(dLu{i}*dLu{j}) ...
            - (Xu*Xv + Yu*Yv)*(dLu{i}*dLv{j} + dLv{i}*dLu{j}) ...
            + (Xu*Xu + Yu*Yu)*(dLv{i}*dLv{j});
        out{cnt} = num/J;                     % K_ij(p) integrand, Eq. (3)
        out{15 + cnt} = Lam{i}*Lam{j}*J;      % volume part of M_ij(p), Eq. (3)
    end
end
for i = 1:5
    out{30 + i} = Lam{i}*J;                   % b_i(p), rank-one mass correction
end
end
