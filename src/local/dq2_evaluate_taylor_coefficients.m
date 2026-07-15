function [CK, CM, Cb] = dq2_evaluate_taylor_coefficients(TC, ea, eb, ec, ed)
% Evaluate all Taylor coefficient matrices at (ea..ed) of any class
% (intval or gradient over intval). Returns cells CK{j+1} 5x5 etc., j=0..9.
% These coefficients are used to enclose the Taylor quotients Delta_t K(e)
% and Delta_t M(e) in Eqs. (29)-(30) and Algorithm 1, Step 1, Eq. (54).
%
% PDF notation:
%   K_j(e), M_j(e), b_j(e) are the symbolic Taylor coefficients used in
%   K(te), M(te), Delta_t K(e), and Delta_t M(e), Eqs. (3), (29)-(30).
%
% Evaluation:
%   If ea, eb, ec, ed are INTLAB intervals, this returns interval enclosures.
%   If they are INTLAB Hessian variables, this returns value/gradient/Hessian
%   data for centered-form evaluation in dq2_algorithm1_box.m.
M = size(TC.MEXP, 1);
% power tables
dmax = max(TC.MEXP(:));
pa = cell(dmax+1,1); pb = pa; pc = pa; pd = pa;
pa{1} = ea*0 + 1; pb{1} = eb*0 + 1; pc{1} = ec*0 + 1; pd{1} = ed*0 + 1;
for k = 1:dmax
    pa{k+1} = pa{k}*ea; pb{k+1} = pb{k}*eb; pc{k+1} = pc{k}*ec; pd{k+1} = pd{k}*ed;
end
mv = cell(M,1);
for r = 1:M
    mv{r} = pa{TC.MEXP(r,1)+1}*pb{TC.MEXP(r,2)+1}*pc{TC.MEXP(r,3)+1}*pd{TC.MEXP(r,4)+1};
end
% stack into a column vector of the right class
monvals = vertcat(mv{:});
CK = cell(10,1); CM = cell(10,1); Cb = cell(10,1);
for j = 0:9
    vK = TC.CK{j+1}*monvals;       % K_j(e), symbolic monomials -> interval/Hessian
    vM = TC.CM{j+1}*monvals;       % M_j(e), symbolic monomials -> interval/Hessian
    vb = TC.Cb{j+1}*monvals;       % b_j(e), symbolic monomials -> interval/Hessian
    K5 = expand5(vK, TC.pr); M5 = expand5(vM, TC.pr);
    CK{j+1} = K5; CM{j+1} = M5; Cb{j+1} = vb;
end
end

function A = expand5(v, pr)
% symmetric 5x5 from 15-vector, preserving class of v
A = [v(1) v(2) v(3) v(4) v(5);
     v(2) v(6) v(7) v(8) v(9);
     v(3) v(7) v(10) v(11) v(12);
     v(4) v(8) v(11) v(13) v(14);
     v(5) v(9) v(12) v(14) v(15)];
end
