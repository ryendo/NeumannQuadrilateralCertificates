function blocks = qn_single_box_blocks(CK,CM,Cb,EK,EM,Eb,ad2,t,FR)
% Raw Taylor-quotient blocks for the local single-box certificate.
%
% This stage deliberately forms no inverse. The caller can therefore
% certify M(te)>0 and the perpendicular block B0 at nu=0 before constructing
% the Schur data used for Fhat_t.

PI2=FR.PI2; V5=FR.V5; W5=FR.W5; D0=FR.D0; I3=FR.I3;

RK=CK{3}; RMv=CM{3}; Rb=Cb{3}; tp=t;
for j=3:9
    RK=RK+tp*CK{j+1}; RMv=RMv+tp*CM{j+1}; Rb=Rb+tp*Cb{j+1};
    tp=tp*t;
end
RK=RK+tp*EK; RMv=RMv+tp*EM; Rb=Rb+tp*Eb;
q_te=1-sqr(t)*ad2; % q(t*e)=|Q_{t e}| in (25)
wv=Cb{2}+t*Rb;
RM=RMv-(wv*wv')/q_te;
Xt=V5'*(RK-PI2*RM)*V5;
DtM=CM{2}+t*RM; DtK=CK{2}+t*RK;
Th=DtK-PI2*DtM;
Yd=V5'*DtM*V5;
Nm=V5'*DtM*W5;
MW=I3+t*(W5'*DtM*W5);
C0=V5'*Th*W5;
B0=D0+t*(W5'*Th*W5);
X0=V5'*(CK{2}-PI2*CM{2})*V5;
al=(X0(1,1)-X0(2,2))/2;
be=(X0(1,2)+X0(2,1))/2;

blocks=struct('Xt',Xt,'DtM',DtM,'DtK',DtK,'Yd',Yd,'Nm',Nm, ...
    'MW',MW,'C0',C0,'B0',B0,'X0',X0,'al',al,'be',be);
end
