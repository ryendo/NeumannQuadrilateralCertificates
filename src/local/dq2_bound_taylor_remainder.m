function [BK, BM, Bb] = dq2_bound_taylor_remainder(Ebox, tb, Ncell, ord)
% Per-box remainder bounds: |order-`ord` Taylor coefficient of the integrals|
% for tau in [0, tb], e in Ebox. Returns sup-magnitude matrices.
% These bounds make the Taylor quotients rigorous on boxes touching t = 0;
% they implement the limiting definitions in Eqs. (29)-(30) and the interval
% enclosures in Algorithm 1, Step 1, Eq. (54).
%
% PDF notation:
%   R_K(t,e), R_M(t,e), R_b(t,e) are implementation remainders for
%   Delta_t K(e), Delta_t M(e), Eqs. (29)-(30).
%
% Evaluation:
%   This is direct INTLAB Taylor arithmetic over tau in [0,tb], e in Ebox,
%   and integration subboxes (u,v).  The output is a radius-only interval
%   bound used as midrad(0, B*) in dq2_algorithm1_box.m.
if nargin < 3, Ncell = 2; end
if nargin < 4, ord = 10; end
PI = intval('pi');
h = intval(1)/intval(Ncell);
edges = intval('-0.5') + h*intval(0:Ncell);
% Accumulate the cellwise upper bounds in interval arithmetic.  Accumulating
% the extracted endpoints in ordinary doubles would not guarantee that the
% sum is rounded upward.
BK = intval(zeros(5)); BM = intval(zeros(5)); Bb = intval(zeros(5,1));
pr = zeros(15,2); c = 0;
for i = 1:5
    for j = i:5
        c = c + 1; pr(c,:) = [i j];
    end
end
tau = taylorinit(infsup(0, sup(intval(tb))), ord); % Taylor variable for p = tau*e
aa = tau*Ebox(1); bb = tau*Ebox(2); cc = tau*Ebox(3); dd = tau*Ebox(4); % p = tau*e
for iu = 1:Ncell
    for iv = 1:Ncell
        uu = infsup(inf(edges(iu)), sup(edges(iu+1)));
        vv = infsup(inf(edges(iv)), sup(edges(iv+1)));
        out = dq2_quadrilateral_kernels(uu, vv, aa, bb, cc, dd, PI);
        for kk = 1:35
            g = out{kk};
            co = g{ord};                      % order-ord Taylor coefficient
            w = sqr(h)*intval(mag(co));       % nonnegative interval upper bound for the cell integral
            if kk <= 15
                BK(pr(kk,1), pr(kk,2)) = BK(pr(kk,1), pr(kk,2)) + w;
            elseif kk <= 30
                BM(pr(kk-15,1), pr(kk-15,2)) = BM(pr(kk-15,1), pr(kk-15,2)) + w;
            else
                Bb(kk-30) = Bb(kk-30) + w;
            end
        end
    end
end
BK = BK + triu(BK,1)'; BM = BM + triu(BM,1)';
% The callers need radius matrices.  Extract only after the complete sum has
% been formed with outward-rounded interval additions.
BK = sup(BK); BM = sup(BM); Bb = sup(Bb);
end
