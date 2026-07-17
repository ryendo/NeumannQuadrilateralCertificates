function [epsK, epsM, info, epsGradK, epsGradM] = qn_gl_pad(p,N)
% Certified Bernstein-ellipse pad for the spatial tensor GL rule.
%
% After exact metric cancellation, both K and M integrands are entire in the
% reference variables: K has the form (q_i'*q_j)*J and M has the form
% phi_i*phi_j*J.  Consequently the pad uses an entire-function ellipse and
% has no Jacobian lower bound or pole-distance assumption.

if nargin<2, N=20; end
PI = intval('pi'); twoPI = intval('2')*PI;
A = sup(abs(p(1))); B = sup(abs(p(2)));
C = sup(abs(p(3))); D = sup(abs(p(4)));
Xu_ub = intval('1') + intval(A) + intval(B);
Yv_ub = intval('1') + intval(A) + intval(C);
Xv_ub = intval(D) + intval(B);
Yu_ub = intval(D) + intval(C);
freq_u = twoPI*(Xu_ub+Yu_ub);
freq_v = twoPI*(Xv_ub+Yv_ub);
[euK,eguK,euM,eguM] = one_direction(freq_u,A,B,C,D,N);
[evK,egvK,evM,egvM] = one_direction(freq_v,A,B,C,D,N);
epsK = intval('2')*(euK+evK);
epsM = intval('2')*(euM+evM);
epsGradK = intval('2')*(eguK+egvK);
epsGradM = intval('2')*(eguM+egvM);
info = struct('representation','physical_gradient_entire', ...
    'freq_u',freq_u,'freq_v',freq_v, ...
    'epsK',epsK,'epsM',epsM, ...
    'epsGradK',epsGradK,'epsGradM',epsGradM);
end

function [eK,egK,eM,egM] = one_direction(freq,A,B,C,D,N)
one = intval('1'); a = intval('1.5');
% A fixed ellipse is deliberately used: any a>1 is valid for these entire
% integrands, and a point interval avoids a parameter-dependent ellipse.
b = sqrt(a^2-one); rho = a+b;

% On the ellipse in one reference coordinate, use |u|,|v| <= a/2 for a
% symmetric conservative bound.  J is affine in u,v.  The four formulas
% below bound all parameter derivatives of J on the same ellipse.
shape = a^2/intval('2');
Jell = one+A^2+D^2 + a*((B+C)*(one+A)+B*D+C*D);
dJa = intval('2')*A+a*(B+C);
dJb = a*(one+A+D);
dJc = dJb;
dJd = intval('2')*D+a*(B+C);
dJ = max(max(dJa,dJb),max(dJc,dJd));

% |q_i'*q_j| <= 8*pi^2 times the combined trig envelope.  A parameter
% derivative of that product is bounded by 64*pi^3*shape; the constants
% cover both physical components and both differentiated factors.
coefK = intval('8')*intval('pi')^2*Jell;
coefgK = intval('64')*intval('pi')^3*shape*Jell ...
    + intval('8')*intval('pi')^2*dJ;
% The same bound covers raw mass entries and single-mode means.
coefM = Jell;
coefgM = intval('4')*intval('pi')*shape*Jell+dJ;
trigEnvelope=cosh(freq*b);
Cgl = intval('64')/intval('15');
factor=Cgl/(rho^(2*N))/(one-one/rho^2);
eK=factor*trigEnvelope*coefK;
egK=factor*trigEnvelope*coefgK;
eM=factor*trigEnvelope*coefM;
egM=factor*trigEnvelope*coefgM;
end
