function [K,Mraw,means] = qn_assemble_ad(p)
% GL assembly for either intval or INTLAB gradient parameter inputs.
%
% The trial functions are restrictions of fixed ambient trigonometric
% functions.  Therefore, if q_i is their physical (X,Y)-gradient, then
% grad_{u,v}(psi_i o Phi)=DPhi'*q_i and the metric factors cancel exactly:
%
%   (grad_ref psi_i)'*(DPhi\DPhi'\grad_ref psi_j)*J = (q_i'*q_j)*J.
%
% Assemble this regular form directly.  In particular, do not form the
% removable quotient 0/0 at a triangular boundary where J vanishes at one
% corner of the reference square.

[x,w]=qn_gauss_legendre_20(); PI=intval('pi'); half=intval('0.5'); two=intval('2');
modes=[1 0;0 1;1 1;2 0;0 2]; z=p(1)*intval('0');
K=cell(5); Mraw=cell(5); means=cell(5,1);
for i=1:5
    means{i}=z;
    for j=1:5, K{i,j}=z; Mraw{i,j}=z; end
end
a=p(1); b=p(2); c=p(3); d=p(4);
for iu=1:20
    u=x(iu);
    for iv=1:20
        v=x(iv); W=w(iu)*w(iv);
        X=u-a*u-d*v-two*b*u*v; Y=v-d*u+a*v+two*c*u*v;
        Xu=intval('1')-a-two*b*v; Xv=-(d+two*b*u);
        Yu=-(d-two*c*v); Yv=intval('1')+a+two*c*u;
        J=Xu*Yv-Xv*Yu;
        phi=cell(5,1); psiX=cell(5,1); psiY=cell(5,1);
        for q=1:5
            mm=intval(modes(q,1)); nn=intval(modes(q,2));
            cx=cos(mm*PI*(X+half)); sx=sin(mm*PI*(X+half));
            cy=cos(nn*PI*(Y+half)); sy=sin(nn*PI*(Y+half));
            phi{q}=cx*cy;
            psiX{q}=-mm*PI*sx*cy;
            psiY{q}=-nn*PI*cx*sy;
        end
        for i=1:5
            means{i}=means{i}+phi{i}*J*W;
            for j=i:5
                K{i,j}=K{i,j}+(psiX{i}*psiX{j}+psiY{i}*psiY{j})*J*W;
                Mraw{i,j}=Mraw{i,j}+phi{i}*phi{j}*J*W;
            end
        end
    end
end
for i=1:5, for j=i+1:5
    K{j,i}=K{i,j}; Mraw{j,i}=Mraw{i,j};
end, end
end
