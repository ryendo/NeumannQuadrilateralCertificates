function [K,Mraw,means] = qn_assemble_interval_center(p)
% Vectorized certified assembly at a point interval in parameter space.
%
% The exactly cancelled physical-gradient form remains regular on triangle
% faces.

[U,V,W]=qn_gl_tensor20();
PI=intval('pi'); half=intval('0.5'); two=intval('2'); one=intval('1');
modes=[1 0;0 1;1 1;2 0;0 2]; nnode=size(U,1);
a=p(1); b=p(2); c=p(3); d=p(4);

UV=U.*V;
X=U-a.*U-d.*V-two.*b.*UV;
Y=V-d.*U+a.*V+two.*c.*UV;
Xu=one-a-two.*b.*V; Xv=-(d+two.*b.*U);
Yu=-(d-two.*c.*V); Yv=one+a+two.*c.*U;
J=Xu.*Yv-Xv.*Yu;

Phi=infsup(zeros(nnode,5),zeros(nnode,5));
PsiX=infsup(zeros(nnode,5),zeros(nnode,5));
PsiY=infsup(zeros(nnode,5),zeros(nnode,5));
for q=1:5
    mx=intval(modes(q,1))*PI; ny=intval(modes(q,2))*PI;
    cx=cos(mx.*(X+half)); sx=sin(mx.*(X+half));
    cy=cos(ny.*(Y+half)); sy=sin(ny.*(Y+half));
    Phi(:,q)=cx.*cy;
    PsiX(:,q)=-mx.*sx.*cy;
    PsiY(:,q)=-ny.*cx.*sy;
end

K=cell(5); Mraw=cell(5); means=cell(5,1);
JW=J.*W;
for i=1:5
    means{i}=sum(Phi(:,i).*JW);
    for j=i:5
        K{i,j}=sum((PsiX(:,i).*PsiX(:,j)+PsiY(:,i).*PsiY(:,j)).*JW);
        Mraw{i,j}=sum(Phi(:,i).*Phi(:,j).*JW);
    end
end
for i=1:5
    for j=i+1:5
        K{j,i}=K{i,j}; Mraw{j,i}=Mraw{i,j};
    end
end
end
