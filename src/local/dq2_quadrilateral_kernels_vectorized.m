function out=dq2_quadrilateral_kernels_vectorized(u,v,aa,bb,cc,dd,PI)
% Certified quadrilateral kernels for vectorized Taylor-remainder cells.
%
% u and v are equally sized intval arrays.  Parameters may be INTLAB Taylor
% objects.  Elementwise operators evaluate the same exact integrands on all
% integration cells at once.

modes=[1 0;0 1;1 1;2 0;0 2];
X=u-aa.*u-dd.*v-intval('2').*bb.*u.*v;
Y=v-dd.*u+aa.*v+intval('2').*cc.*u.*v;
Xu=intval('1')-aa-intval('2').*bb.*v;
Xv=-(dd+intval('2').*bb.*u);
Yu=-(dd-intval('2').*cc.*v);
Yv=intval('1')+aa+intval('2').*cc.*u;
J=Xu.*Yv-Xv.*Yu;

Lam=cell(5,1); dLu=cell(5,1); dLv=cell(5,1);
for i=1:5
    m=intval(modes(i,1)); n=intval(modes(i,2));
    cx=cos(m.*PI.*(X+intval('0.5'))); sx=sin(m.*PI.*(X+intval('0.5')));
    cy=cos(n.*PI.*(Y+intval('0.5'))); sy=sin(n.*PI.*(Y+intval('0.5')));
    Lam{i}=cx.*cy;
    Lx=-m.*PI.*sx.*cy; Ly=-n.*PI.*cx.*sy;
    dLu{i}=Lx.*Xu+Ly.*Yu;
    dLv{i}=Lx.*Xv+Ly.*Yv;
end

out=cell(35,1); cnt=0;
for i=1:5
    for j=i:5
        cnt=cnt+1;
        num=(Xv.*Xv+Yv.*Yv).*dLu{i}.*dLu{j} ...
            -(Xu.*Xv+Yu.*Yv).*(dLu{i}.*dLv{j}+dLv{i}.*dLu{j}) ...
            +(Xu.*Xu+Yu.*Yu).*dLv{i}.*dLv{j};
        out{cnt}=num./J;
        out{15+cnt}=Lam{i}.*Lam{j}.*J;
    end
end
for i=1:5, out{30+i}=Lam{i}.*J; end
end
