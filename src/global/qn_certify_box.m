function [verdict,info,gap,split_dim] = qn_certify_box(center,half_widths)
% Certified per-box Rayleigh bound from Proposition p:box-bound.
% Floating-point work chooses only the frame and subdivision coordinate.

C=qn_global_constants(); center=center(:); half_widths=half_widths(:);
[Kc,Mc]=qn_km_float(center);
[V,D]=eig((Kc+Kc')/2,(Mc+Mc')/2);
[evals,idx]=sort(real(diag(D))); V=real(V(:,idx));
gap=evals(2)-evals(1); v=V(:,1);

% Slack-driven coordinate, with longest-side fallback. This is non-certified.
slack=zeros(4,1);
for k=1:4
    h=1e-6*max(1,abs(center(k))); pp=center; pm=center;
    pp(k)=pp(k)+h; pm(k)=pm(k)-h;
    [Kp,Mp]=qn_km_float(pp); [Km,Mm]=qn_km_float(pm);
    dK=(Kp-Km)/(2*h); dM=(Mp-Mm)/(2*h);
    slack(k)=(abs(v'*dK*v)+abs(v'*dM*v))*half_widths(k);
end
[mx,split_dim]=max(slack);
[longest, longest_dim]=max(half_widths);
if ~(isfinite(mx) && mx>0) || half_widths(split_dim)<longest/2
    split_dim=longest_dim;
end

try
    [K,M,area]=qn_km_enclosure(center,half_widths);
catch ME
    if startsWith(ME.identifier,'qn:')
        verdict='fail'; info=struct('reason',ME.identifier); return;
    end
    rethrow(ME);
end
if ~qn_interval_ldl_pd(M)
    verdict='fail'; info=struct('reason','mass_not_pd'); return;
end
qh=intval(sup(area));

if gap<C.gap_threshold
    [ok,info]=try_two_vector(V(:,1:2),K,M,qh,C.pi2);
    if ok, verdict='cert'; return; end
end
N=quad_form(v,K); Den=quad_form(v,M);
if inf(Den)>0
    ub=intval(sup(N))/intval(inf(Den));
    if sup(qh*ub)<inf(C.pi2)
        verdict='cert'; info=struct('route','1v','margin',inf(C.pi2-qh*ub)); return;
    end
end
if gap>=C.gap_threshold
    [ok,info]=try_two_vector(V(:,1:2),K,M,qh,C.pi2);
    if ok, verdict='cert'; return; end
end
verdict='fail'; info=struct('reason','rayleigh_bound');
end

function [ok,out]=try_two_vector(V,K,M,qh,pi2)
ub=two_vector_upper_bound(V,K,M);
ok=~isempty(ub) && sup(qh*ub)<inf(pi2);
if ok
    out=struct('route','2v','margin',inf(pi2-qh*ub));
else
    out=struct();
end
end

function q=quad_form(v,A)
vi=intval(zeros(numel(v),1));
for k=1:numel(v), vi(k)=intval(sprintf('%.17g',v(k))); end
q=vi'*A*vi;
end

function ub=two_vector_upper_bound(V,K,M)
Vi=intval(zeros(size(V)));
for i=1:size(V,1), for j=1:2
    Vi(i,j)=intval(sprintf('%.17g',V(i,j)));
end, end
A=Vi'*K*Vi; C=Vi'*M*Vi;
detC=C(1,1)*C(2,2)-C(1,2)^2;
if inf(C(1,1))<=0 || inf(detC)<=0, ub=[]; return; end
beta=A(1,1)*C(2,2)+A(2,2)*C(1,1)-intval('2')*A(1,2)*C(1,2);
gamma=A(1,1)*A(2,2)-A(1,2)^2;
disc=beta^2-intval('4')*detC*gamma;
if inf(disc)>=0
    root=(beta-sqrt(disc))/(intval('2')*detC);
    ub=intval(sup(root)); return;
end
best=inf;
for k=0:89
    th=pi*k/89; ct=intval(sprintf('%.17g',cos(th))); st=intval(sprintf('%.17g',sin(th)));
    num=ct^2*A(1,1)+intval('2')*ct*st*A(1,2)+st^2*A(2,2);
    den=ct^2*C(1,1)+intval('2')*ct*st*C(1,2)+st^2*C(2,2);
    if inf(den)>0
        cand=intval(sup(num))/intval(inf(den));
        best=min(best,sup(cand));
    end
end
if isfinite(best), ub=intval(best); else, ub=[]; end
end
