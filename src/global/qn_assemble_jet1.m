function [K,Mraw,means] = qn_assemble_jet1(p)
% Lightweight interval forward AD: each jet stores value and 4-gradient.
% Assemble the exactly cancelled physical-gradient form (q_i'*q_j)*J;
% see qn_assemble_ad.  This form remains regular when J=0 at a corner.

[x,w]=qn_gauss_legendre_20(); PI=intval('pi'); half=intval('0.5'); two=intval('2');
modes=[1 0;0 1;1 1;2 0;0 2]; z=jet(intval('0'),intval(zeros(4,1)));
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
        X=jet(u-a*u-d*v-two*b*u*v,[-u;-two*u*v;intval('0');-v]);
        Y=jet(v-d*u+a*v+two*c*u*v,[v;intval('0');two*u*v;-u]);
        Xu=jet(intval('1')-a-two*b*v,[-intval('1');-two*v;intval('0');intval('0')]);
        Xv=jet(-(d+two*b*u),[intval('0');-two*u;intval('0');-intval('1')]);
        Yu=jet(-(d-two*c*v),[intval('0');intval('0');two*v;-intval('1')]);
        Yv=jet(intval('1')+a+two*c*u,[intval('1');intval('0');two*u;intval('0')]);
        J=jsub(jmul(Xu,Yv),jmul(Xv,Yu));
        phi=cell(5,1); psiX=cell(5,1); psiY=cell(5,1);
        for q=1:5
            mm=intval(modes(q,1)); nn=intval(modes(q,2));
            ax=jscale(jaddc(X,half),mm*PI); ay=jscale(jaddc(Y,half),nn*PI);
            cx=jcos(ax); sx=jsin(ax); cy=jcos(ay); sy=jsin(ay);
            phi{q}=jmul(cx,cy);
            psiX{q}=jscale(jmul(sx,cy),-mm*PI);
            psiY{q}=jscale(jmul(cx,sy),-nn*PI);
        end
        for i=1:5
            means{i}=jadd(means{i},jscale(jmul(phi{i},J),W));
            for j=i:5
                qdot=jadd(jmul(psiX{i},psiX{j}),jmul(psiY{i},psiY{j}));
                K{i,j}=jadd(K{i,j},jscale(jmul(qdot,J),W));
                Mraw{i,j}=jadd(Mraw{i,j},jscale(jmul(jmul(phi{i},phi{j}),J),W));
            end
        end
    end
end
for i=1:5, for j=i+1:5
    K{j,i}=K{i,j}; Mraw{j,i}=Mraw{i,j};
end, end
end

function z=jet(v,g), z=struct('v',v,'g',g); end
function z=jadd(x,y), z=jet(x.v+y.v,x.g+y.g); end
function z=jsub(x,y), z=jet(x.v-y.v,x.g-y.g); end
function z=jaddc(x,c), z=jet(x.v+c,x.g); end
function z=jscale(x,c), z=jet(x.v*c,x.g*c); end
function z=jmul(x,y), z=jet(x.v*y.v,x.g*y.v+x.v*y.g); end
function z=jsin(x), z=jet(sin(x.v),cos(x.v)*x.g); end
function z=jcos(x), z=jet(cos(x.v),-sin(x.v)*x.g); end
