function [K,M,q] = qn_km_float(p)
% Floating-point center assembly used only to choose the coordinate in (47).
% The certified box decision is made independently by qn_km_enclosure and
% qn_certify_box.  The formulas are the fixed-reference-square integrals of
% Section 4.3, with q(p)=|Q_p| as in equation (2).

persistent U V W
if isempty(U)
    [xi,wi]=qn_gauss_legendre_20(); x=mid(xi); w=mid(wi);
    Ug=repmat(x,1,20); Vg=repmat(x.',20,1);
    Wg=repmat(w,1,20).*repmat(w.',20,1);
    U=Ug(:); V=Vg(:); W=Wg(:);
end
modes=[1 0; 0 1; 1 1; 2 0; 0 2];
a=p(1); b=p(2); c=p(3); d=p(4);
UV=U.*V;
X=U-a.*U-d.*V-2*b.*UV;
Y=V-d.*U+a.*V+2*c.*UV;
Xu=1-a-2*b.*V; Xv=-(d+2*b.*U);
Yu=-(d-2*c.*V); Yv=1+a+2*c.*U;
J=Xu.*Yv-Xv.*Yu; JW=J.*W; q=sum(JW);

Phi=zeros(400,5); PsiX=zeros(400,5); PsiY=zeros(400,5);
for k=1:5
    m=modes(k,1); n=modes(k,2);
    cx=cos(m*pi*(X+0.5)); sx=sin(m*pi*(X+0.5));
    cy=cos(n*pi*(Y+0.5)); sy=sin(n*pi*(Y+0.5));
    Phi(:,k)=cx.*cy;
    PsiX(:,k)=-m*pi.*sx.*cy;
    PsiY(:,k)=-n*pi.*cx.*sy;
end

K=zeros(5); R=zeros(5); means=zeros(5,1);
for i=1:5
    means(i)=sum(Phi(:,i).*JW);
    for j=i:5
        K(i,j)=sum((PsiX(:,i).*PsiX(:,j)+ ...
            PsiY(:,i).*PsiY(:,j)).*JW);
        R(i,j)=sum(Phi(:,i).*Phi(:,j).*JW);
    end
end
K=K+triu(K,1)'; R=R+triu(R,1)';
M=R-(means*means')/q;
end
