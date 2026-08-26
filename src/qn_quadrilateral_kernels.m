function out=qn_quadrilateral_kernels(u,v,a,b,c,d)
% Exact stiffness, raw-mass, and mean integrands for the five trial modes.

PI=intval('pi'); two=intval('2'); half=intval('0.5');
modes=[1 0;0 1;1 1;2 0;0 2];
X=u-a.*u-d.*v-two.*b.*u.*v;
Y=v-d.*u+a.*v+two.*c.*u.*v;
Xu=intval('1')-a-two.*b.*v; Xv=-(d+two.*b.*u);
Yu=-(d-two.*c.*v); Yv=intval('1')+a+two.*c.*u;
J=Xu.*Yv-Xv.*Yu;

phi=cell(5,1); psiX=cell(5,1); psiY=cell(5,1);
for i=1:5
    mx=intval(modes(i,1))*PI; ny=intval(modes(i,2))*PI;
    cx=cos(mx.*(X+half)); sx=sin(mx.*(X+half));
    cy=cos(ny.*(Y+half)); sy=sin(ny.*(Y+half));
    phi{i}=cx.*cy;
    psiX{i}=-mx.*sx.*cy;
    psiY{i}=-ny.*cx.*sy;
end

out=cell(35,1); k=0;
for i=1:5
    for j=i:5
        k=k+1;
        out{k}=(psiX{i}.*psiX{j}+psiY{i}.*psiY{j}).*J;
        out{15+k}=phi{i}.*phi{j}.*J;
    end
end
for i=1:5, out{30+i}=phi{i}.*J; end
end
