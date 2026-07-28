function [BK,BM,Bb]=dq2_bound_taylor_remainder_vectorized(Ebox,tb,Ncell,ord)
% Vectorized certified Taylor-remainder bound over all integration cells.

if nargin<3, Ncell=2; end
if nargin<4, ord=10; end
PI=intval('pi'); h=intval('1')/intval(Ncell);
edges=intval('-0.5')+h*intval(0:Ncell);
ncell=Ncell^2; ulo=zeros(ncell,1); uhi=ulo; vlo=ulo; vhi=ulo; r=0;
for iu=1:Ncell
    for iv=1:Ncell
        r=r+1;
        ulo(r)=inf(edges(iu)); uhi(r)=sup(edges(iu+1));
        vlo(r)=inf(edges(iv)); vhi(r)=sup(edges(iv+1));
    end
end
u=infsup(ulo,uhi); v=infsup(vlo,vhi);

tau=taylorinit(infsup(0,sup(intval(tb))),ord);
aa=tau*Ebox(1); bb=tau*Ebox(2); cc=tau*Ebox(3); dd=tau*Ebox(4);
out=dq2_quadrilateral_kernels_vectorized(u,v,aa,bb,cc,dd,PI);

BK=intval(zeros(5)); BM=intval(zeros(5)); Bb=intval(zeros(5,1));
pr=zeros(15,2); cpair=0;
for i=1:5
    for j=i:5, cpair=cpair+1; pr(cpair,:)=[i j]; end
end
h2=sqr(h);
for kk=1:35
    co=out{kk}{ord};
    radius=h2*sum(intval(mag(co)));
    if kk<=15
        BK(pr(kk,1),pr(kk,2))=radius;
    elseif kk<=30
        BM(pr(kk-15,1),pr(kk-15,2))=radius;
    else
        Bb(kk-30)=radius;
    end
end
BK=BK+triu(BK,1)'; BM=BM+triu(BM,1)';
BK=sup(BK); BM=sup(BM); Bb=sup(Bb);
end
