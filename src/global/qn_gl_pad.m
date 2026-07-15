function [epsK, epsM, info] = qn_gl_pad(p,N)
% Certified Bernstein-ellipse pad for the 12x12 spatial GL rule.
% Mechanical INTLAB translation of _kernel/gl_pad_v2.py (value slot).

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
g_ub = max((Xu_ub+Yu_ub)^2, (Xv_ub+Yv_ub)^2);
kcoef = twoPI^2*g_ub;
Jhi = Xu_ub*Yv_ub + Xv_ub*Yu_ub;

half = intval('0.5');
J00 = jac(-half,-half,p); J10 = jac(half,-half,p);
J01 = jac(-half,half,p);  J11 = jac(half,half,p);
minJ = intval(min([inf(J00),inf(J10),inf(J01),inf(J11)]));
if inf(minJ) <= 0, error('qn:Jacobian','Jacobian is not certifiably positive.'); end
dJdu = intval(sup(abs(J10-J00)));
dJdv = intval(sup(abs(J01-J00)));

euK = one_direction(freq_u,dJdu,minJ,kcoef,Jhi,true,N);
evK = one_direction(freq_v,dJdv,minJ,kcoef,Jhi,true,N);
euM = one_direction(freq_u,dJdu,minJ,kcoef,Jhi,false,N);
evM = one_direction(freq_v,dJdv,minJ,kcoef,Jhi,false,N);
epsK = intval('2')*(euK+evK);
epsM = intval('2')*(euM+evM);
info = struct('minJ',minJ,'freq_u',freq_u,'freq_v',freq_v, ...
    'epsK',epsK,'epsM',epsM);
end

function e = one_direction(freq,dJ,minJ,kcoef,Jhi,isK,N)
one = intval('1'); cap = intval('1.5');
if sup(dJ) <= sup(intval('1e-9')*minJ)
    xstar = intval('1e6');
else
    xstar = one + intval('2')*minJ/dJ;
end
if inf(xstar) <= 1, error('qn:GLPole','GL pole reach is not greater than one.'); end
a = one + (xstar-one)/intval('2');
a = min(a,cap);
b = sqrt(a^2-one); rho = a+b;
if isK
    poly = kcoef*(one+a+a^2+a^3+a^4);
    num = cosh(freq*b)*poly;
    Jell = minJ-dJ*a;
    if inf(Jell) <= 0, Jell = minJ/intval('2'); end
    Mrho = num/Jell;
else
    poly = Jhi*(one+a);
    Mrho = cosh(freq*b)*poly;
end
Cgl = intval('64')/intval('15');
e = Cgl*Mrho/(rho^(2*N))/(one-one/rho^2);
end

function J = jac(u,v,p)
a=p(1); b=p(2); c=p(3); d=p(4);
Xu=intval('1')-a-intval('2')*b*v;
Xv=-(d+intval('2')*b*u);
Yu=-(d-intval('2')*c*v);
Yv=intval('1')+a+intval('2')*c*u;
J=Xu*Yv-Xv*Yu;
end
