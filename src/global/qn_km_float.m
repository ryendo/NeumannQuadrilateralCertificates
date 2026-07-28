function [K, M, area] = qn_km_float(p)
% Non-certified center assembly. Used only to choose a test frame and split axis.
% Use the same exactly cancelled physical-gradient stiffness form as the
% certified assembly, so center evaluation is also defined on triangle faces.

persistent U V W
if isempty(U)
    [xi,wi]=qn_gauss_legendre_20(); x=mid(xi); w=mid(wi);
    Ug=repmat(x,1,20); Vg=repmat(x.',20,1);
    Wg=repmat(w,1,20).*repmat(w.',20,1);
    U=Ug(:); V=Vg(:); W=Wg(:);
end
modes = [1 0; 0 1; 1 1; 2 0; 0 2];
a = p(1); b = p(2); c = p(3); d = p(4);
UV=U.*V;
X=U-a.*U-d.*V-2*b.*UV;
Y=V-d.*U+a.*V+2*c.*UV;
Xu=1-a-2*b.*V; Xv=-(d+2*b.*U);
Yu=-(d-2*c.*V); Yv=1+a+2*c.*U;
J=Xu.*Yv-Xv.*Yu; JW=J.*W; area=sum(JW);

Phi=zeros(400,5); PsiX=zeros(400,5); PsiY=zeros(400,5);
for q=1:5
    mm=modes(q,1); nn=modes(q,2);
    cx=cos(mm*pi*(X+0.5)); sx=sin(mm*pi*(X+0.5));
    cy=cos(nn*pi*(Y+0.5)); sy=sin(nn*pi*(Y+0.5));
    Phi(:,q)=cx.*cy;
    PsiX(:,q)=-mm*pi.*sx.*cy;
    PsiY(:,q)=-nn*pi.*cx.*sy;
end

K=zeros(5); Mraw=zeros(5); means=zeros(5,1);
for i=1:5
    means(i)=sum(Phi(:,i).*JW);
    for j=i:5
        K(i,j)=sum((PsiX(:,i).*PsiX(:,j)+PsiY(:,i).*PsiY(:,j)).*JW);
        Mraw(i,j)=sum(Phi(:,i).*Phi(:,j).*JW);
    end
end
K = K + triu(K,1)'; Mraw = Mraw + triu(Mraw,1)';
M = Mraw - (means*means')/area;
end
