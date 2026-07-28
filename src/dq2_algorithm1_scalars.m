function [d1_over_t,d2,d0,S_core,schur] = dq2_algorithm1_scalars(CK,CM,Cb,EK,EM,Eb,q,t,FR)
% Schur/Vieta quantities shared by full DQ2 checks and sign-only retests.
% The operation order matches Algorithm 1, Eqs. (29)--(30), (59)--(60), (67).

PI2=FR.PI2; V5=FR.V5; W5=FR.W5; D0=FR.D0; I2=FR.I2; I3=FR.I3;

% Taylor quotient for Delta_t K, Delta_t M, and the mean vector.
RK=CK{3}; RMv=CM{3}; Rb=Cb{3}; tp=t;
for j=3:9
    RK=RK+tp*CK{j+1}; RMv=RMv+tp*CM{j+1}; Rb=Rb+tp*Cb{j+1};
    tp=tp*t;
end
RK=RK+tp*EK; RMv=RMv+tp*EM; Rb=Rb+tp*Eb;
den=1-sqr(t)*q;
wv=Cb{2}+t*Rb;
RM=RMv-(wv*wv')/den;
Xt=V5'*(RK-PI2*RM)*V5;
DtM=CM{2}+t*RM; DtK=CK{2}+t*RK;
Th=DtK-PI2*DtM;
Yd=V5'*DtM*V5;
Nm=V5'*DtM*W5;
MW=I3+t*(W5'*DtM*W5);
C0=V5'*Th*W5;
B0=D0+t*(W5'*Th*W5);
[B0i,ok]=inverse3_if_regular(B0);
if ~ok
    d1_over_t=[]; d2=[]; d0=[]; S_core=[]; schur=[];
    return
end
Z0=symmetrize(C0*B0i*C0');
CB0=C0*B0i;
Zb1=symmetrize(-(Nm*B0i*C0'+C0*B0i*Nm')+CB0*MW*B0i*C0');

% Reduced two-dimensional pencil and its Vieta coefficients.
X0=V5'*(CK{2}-PI2*CM{2})*V5;
al=(X0(1,1)-X0(2,2))/2;
be=(X0(1,2)+X0(2,1))/2;
G=Xt-Z0;
CBm=Yd+t*Zb1;
d2=determinant2(I2+t*CBm);
d1_over_t=(G(1,1)+G(2,2))+al*(CBm(2,2)-CBm(1,1))-2*be*CBm(1,2) ...
    +t*(G(1,1)*CBm(2,2)+G(2,2)*CBm(1,1)-2*G(1,2)*CBm(1,2));
d0=-(class_square(al)+class_square(be)) ...
    +t*(al*(G(2,2)-G(1,1))-2*be*G(1,2)) ...
    +sqr(t)*(G(1,1)*G(2,2)-class_square(G(1,2)));
S_core=-(PI2*d1_over_t+2*d0)/d2 ...
    +2*q*(sqr(PI2)+sqr(t)*(PI2*d1_over_t+d0)/d2);
if nargout>4
    schur=struct('X0',X0,'al',al,'be',be,'Xt',Xt,'DtK',DtK,'DtM',DtM, ...
        'Y',I2+t*Yd,'Nm',Nm,'MW',MW,'C0',C0,'B0',B0,'B0i',B0i, ...
        'Bbar',I2+t*CBm,'Z0',Z0,'G',G);
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
