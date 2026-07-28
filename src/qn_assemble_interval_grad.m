function [K,Mraw,means] = qn_assemble_interval_grad(p)
% Vectorized exact first-parameter derivatives over an interval box.
%
% Each returned entry has fields v and g.  The derivatives below are an exact
% symbolic chain-rule expansion evaluated with INTLAB interval arithmetic.

[U,V,W]=qn_gl_tensor20();
PI=intval('pi'); half=intval('0.5'); two=intval('2'); one=intval('1');
modes=[1 0;0 1;1 1;2 0;0 2]; nnode=size(U,1); npar=4;
a=p(1); b=p(2); c=p(3); d=p(4);
Z=infsup(zeros(nnode,1),zeros(nnode,1));
O=infsup(ones(nnode,1),ones(nnode,1)); UV=U.*V;

X=U-a.*U-d.*V-two.*b.*UV;
Y=V-d.*U+a.*V+two.*c.*UV;
Xu=one-a-two.*b.*V; Xv=-(d+two.*b.*U);
Yu=-(d-two.*c.*V); Yv=one+a+two.*c.*U;

% Exact parameter derivatives, in the order (a,b,c,d).
dX=[-U,-two.*UV,Z,-V];
dY=[V,Z,two.*UV,-U];
dXu=[-O,-two.*V,Z,Z];
dXv=[Z,-two.*U,Z,-O];
dYu=[Z,Z,two.*V,-O];
dYv=[O,Z,two.*U,Z];

J=Xu.*Yv-Xv.*Yu;
dJ=dXu.*repeat_columns(Yv,npar)+repeat_columns(Xu,npar).*dYv ...
    -dXv.*repeat_columns(Yu,npar)-repeat_columns(Xv,npar).*dYu;

phi=cell(5,1); psiX=cell(5,1); psiY=cell(5,1);
dphi=cell(5,1); dpsiX=cell(5,1); dpsiY=cell(5,1);
for q=1:5
    mx=intval(modes(q,1))*PI; ny=intval(modes(q,2))*PI;
    cx=cos(mx.*(X+half)); sx=sin(mx.*(X+half));
    cy=cos(ny.*(Y+half)); sy=sin(ny.*(Y+half));
    cxcy=cx.*cy; sxsy=sx.*sy;
    phi{q}=cxcy;
    psiX{q}=-mx.*sx.*cy;
    psiY{q}=-ny.*cx.*sy;
    dphi{q}=row_scale(dX,-mx.*sx.*cy)+row_scale(dY,-ny.*cx.*sy);
    dpsiX{q}=row_scale(dX,-mx^2.*cxcy)+row_scale(dY,mx.*ny.*sxsy);
    dpsiY{q}=row_scale(dX,mx.*ny.*sxsy)+row_scale(dY,-ny^2.*cxcy);
end

K=cell(5); Mraw=cell(5); means=cell(5,1);
JW=J.*W; W4=repeat_columns(W,npar);
for i=1:5
    meanValue=sum(phi{i}.*JW);
    meanGrad=sum((row_scale(dphi{i},J)+row_scale(dJ,phi{i})).*W4,1).';
    means{i}=make_jet(meanValue,meanGrad);
    for j=i:5
        qdot=psiX{i}.*psiX{j}+psiY{i}.*psiY{j};
        dqdot=row_scale(dpsiX{i},psiX{j})+row_scale(dpsiX{j},psiX{i}) ...
            +row_scale(dpsiY{i},psiY{j})+row_scale(dpsiY{j},psiY{i});
        kValue=sum(qdot.*JW);
        kGrad=sum((row_scale(dqdot,J)+row_scale(dJ,qdot)).*W4,1).';

        phiprod=phi{i}.*phi{j};
        dphiprod=row_scale(dphi{i},phi{j})+row_scale(dphi{j},phi{i});
        mValue=sum(phiprod.*JW);
        mGrad=sum((row_scale(dphiprod,J)+row_scale(dJ,phiprod)).*W4,1).';
        K{i,j}=make_jet(kValue,kGrad);
        Mraw{i,j}=make_jet(mValue,mGrad);
    end
end
for i=1:5
    for j=i+1:5
        K{j,i}=K{i,j}; Mraw{j,i}=Mraw{i,j};
    end
end
end

function Y=row_scale(X,s)
Y=X.*repeat_columns(s,size(X,2));
end

function Y=repeat_columns(x,n)
Y=infsup(repmat(inf(x),1,n),repmat(sup(x),1,n));
end

function z=make_jet(v,g)
z=struct('v',v,'g',g);
end
