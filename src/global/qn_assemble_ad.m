function [K,Mraw,means] = qn_assemble_ad(p)
% GL assembly for either intval or INTLAB gradient parameter inputs.

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
        if isintval(J), Jx=J; else, Jx=J.x; end
        if inf(Jx)<=0, error('qn:Jacobian','AD assembly Jacobian crosses zero.'); end
        gvv=Xv^2+Yv^2; guv=Xu*Xv+Yu*Yv; guu=Xu^2+Yu^2;
        phi=cell(5,1); phiu=cell(5,1); phiv=cell(5,1);
        for q=1:5
            mm=intval(modes(q,1)); nn=intval(modes(q,2));
            cx=cos(mm*PI*(X+half)); sx=sin(mm*PI*(X+half));
            cy=cos(nn*PI*(Y+half)); sy=sin(nn*PI*(Y+half));
            phi{q}=cx*cy;
            psiX=-mm*PI*sx*cy; psiY=-nn*PI*cx*sy;
            phiu{q}=psiX*Xu+psiY*Yu; phiv{q}=psiX*Xv+psiY*Yv;
        end
        for i=1:5
            means{i}=means{i}+phi{i}*J*W;
            for j=i:5
                K{i,j}=K{i,j}+(gvv*phiu{i}*phiu{j} ...
                    -guv*(phiu{i}*phiv{j}+phiv{i}*phiu{j}) ...
                    +guu*phiv{i}*phiv{j})/J*W;
                Mraw{i,j}=Mraw{i,j}+phi{i}*phi{j}*J*W;
            end
        end
    end
end
for i=1:5, for j=i+1:5
    K{j,i}=K{i,j}; Mraw{j,i}=Mraw{i,j};
end, end
end
