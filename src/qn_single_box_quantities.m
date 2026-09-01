function [d1_over_t,d2,d0,S_core,schur] = qn_single_box_quantities(blocks,ad2,t,FR)
% Schur/Vieta quantities shared by the single-box certificate and retests.
% Supplies Algorithm 1 with the affine Schur/Vieta data from (28) and
% (31)-(35), the central parts of the enclosures in (36), and the
% substitution in S from (23). The caller chooses [nu], subdivides it, and
% certifies the remaining conditions for F_t.
% Preconditions: the whole-box caller has certified M(te)>0 and blocks.B0>0.
% This routine is the first stage that forms a Schur inverse.
% ad2=e_a^2+e_d^2.

PI2=FR.PI2; I2=FR.I2;
Xt=blocks.Xt; DtM=blocks.DtM; DtK=blocks.DtK;
Yd=blocks.Yd; Nm=blocks.Nm; MW=blocks.MW;
C0=blocks.C0; B0=blocks.B0;
[B0i,ok]=inverse3_if_regular(B0);
if ~ok
    d1_over_t=[]; d2=[]; d0=[]; S_core=[]; schur=[];
    return
end
Z0=symmetrize(C0*B0i*C0');
CB0=C0*B0i;
Zb1=symmetrize(-(Nm*B0i*C0'+C0*B0i*Nm')+CB0*MW*B0i*C0');

% Reduced two-dimensional pencil and its Vieta coefficients.
X0=blocks.X0; al=blocks.al; be=blocks.be;
G=Xt-Z0;
CBm=Yd+t*Zb1;
minus_dFhat_dnu=I2+t*CBm; % -partial_nu Fhat_t in (35)
d2=determinant2(minus_dFhat_dnu);
d1_over_t=(G(1,1)+G(2,2))+al*(CBm(2,2)-CBm(1,1))-2*be*CBm(1,2) ...
    +t*(G(1,1)*CBm(2,2)+G(2,2)*CBm(1,1)-2*G(1,2)*CBm(1,2));
d0=-(class_square(al)+class_square(be)) ...
    +t*(al*(G(2,2)-G(1,1))-2*be*G(1,2)) ...
    +sqr(t)*(G(1,1)*G(2,2)-class_square(G(1,2)));
S_core=-(PI2*d1_over_t+2*d0)/d2 ...
    +2*ad2*(sqr(PI2)+sqr(t)*(PI2*d1_over_t+d0)/d2);
if nargout>4
    schur=struct('X0',X0,'al',al,'be',be,'Xt',Xt,'DtK',DtK,'DtM',DtM, ...
        'Y',I2+t*Yd,'Nm',Nm,'MW',MW,'C0',C0,'B0',B0,'B0i',B0i, ...
        'minus_dFhat_dnu',minus_dFhat_dnu,'Z0',Z0,'G',G);
end
end

function y=class_square(x)
if isintval(x), y=sqr(x); else, y=x*x; end
end

function S=symmetrize(A)
S=(A+A')/2;
end

function d=determinant2(A)
d=A(1,1)*A(2,2)-A(1,2)*A(2,1);
end

function [Ai,ok]=inverse3_if_regular(A)
D=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) ...
 -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) ...
 +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1));
if isintval(D), Dx=D; else, Dx=D.x; end
ok=inf(Dx)>0 || sup(Dx)<0;
if ~ok, Ai=A; return; end
Ai=[A(2,2)*A(3,3)-A(2,3)*A(3,2), A(1,3)*A(3,2)-A(1,2)*A(3,3), A(1,2)*A(2,3)-A(1,3)*A(2,2);
    A(2,3)*A(3,1)-A(2,1)*A(3,3), A(1,1)*A(3,3)-A(1,3)*A(3,1), A(1,3)*A(2,1)-A(1,1)*A(2,3);
    A(2,1)*A(3,2)-A(2,2)*A(3,1), A(1,2)*A(3,1)-A(1,1)*A(3,2), A(1,1)*A(2,2)-A(1,2)*A(2,1)]/D;
end
