function report = qn_smoke_test(repo_root)
% Quick integration check; requires MATLAB and INTLAB but no full recomputation.

if nargin<1, repo_root=fileparts(fileparts(mfilename('fullpath'))); end
addpath(fullfile(repo_root,'src'));
C=qn_global_constants();
assert(in(3232/(27*pi^6),C.rho_local));
[K0,M0,A0]=qn_km_enclosure(zeros(4,1),zeros(4,1));
PI=intval('pi');
Kexact=diag([PI^2/intval('2'),PI^2/intval('2'),PI^2/intval('2'), ...
    intval('2')*PI^2,intval('2')*PI^2]);
Mexact=diag([intval('0.5'),intval('0.5'),intval('0.25'),intval('0.5'),intval('0.5')]);
assert(all(all(in(mid(Kexact),K0))) && all(all(in(mid(Mexact),M0))) && in(1,A0), ...
    'Square pencil enclosure misses its analytic value.');
center=[0.10;0.05;0.04;0.02]; hw=0.01*ones(4,1);
[K,M,area]=qn_km_enclosure(center,hw); [Kf,Mf,af]=qn_km_float(center);
containsK=all(all(in(Kf,K))); containsM=all(all(in(Mf,M))); containsArea=in(af,area);
assert(containsK && containsM && containsArea,'Certified enclosure misses its center value.');
[verdict,info,gap,sdim]=qn_certify_box(center,hw);
assert(strcmp(verdict,'cert') && strcmp(info.route,'veigs') && ...
    any(info.index_range==1),'Representative box was not certified by veigs.');

% A box crossing a genuine triangle face c_3=0 must be assemblable. The
% physical-gradient stiffness formula has no division by J, even though the
% reference-map Jacobian interval contains zero at the degenerate corner.
a=-0.17123475299635676; b=0.20062825348897717; d=-0.022264254417977088;
c=(1-a*a-d*d-b*(1-d+a))/(1-d-a);
boundaryCenter=[a;b;c;d]; boundaryHw=1e-6*ones(4,1);
pb=qn_interval_box(boundaryCenter,boundaryHw);
c3=intval('1')-pb(1)^2-pb(4)^2-(pb(2)+pb(3))*(intval('1')-pb(4)) ...
    -pb(1)*(pb(2)-pb(3));
assert(inf(c3)<=0 && sup(c3)>=0,'Boundary smoke box does not cross c_3=0.');
[Kb,Mb,Ab]=qn_km_enclosure(boundaryCenter,boundaryHw);
[Kbf,Mbf,Abf]=qn_km_float(boundaryCenter);
boundaryContains=all(all(in(Kbf,Kb))) && all(all(in(Mbf,Mb))) && in(Abf,Ab);
assert(boundaryContains,'Triangle-boundary pencil enclosure misses its center value.');
[boundaryVerdict,boundaryInfo]=qn_certify_box(boundaryCenter,boundaryHw);
assert(strcmp(boundaryVerdict,'cert'), ...
    'Representative triangle-boundary box was not certified.');
assert(strcmp(boundaryInfo.route,'veigs') && any(boundaryInfo.index_range==1), ...
    'Triangle-boundary box was not certified by veigs with index 1.');

local=qn_summarize_local_results(fullfile(repo_root,'results','local'));
assert(local.verified,'Saved local certificate is incomplete.');
report=struct('constants',struct('rho_local',mid(C.rho_local),'rho_seam',mid(C.rho_seam)), ...
    'enclosure',struct('K',containsK,'M',containsM,'area',containsArea), ...
    'triangle_boundary_enclosure',struct('crosses_c3_zero',true, ...
        'contains_center',boundaryContains,'verdict',boundaryVerdict,'info',boundaryInfo), ...
    'sample_box',struct('verdict',verdict,'info',info,'gap',gap,'split_dim',sdim), ...
    'local_saved_results_verified',local.verified);
if isfield(info,'reason'), detail=info.reason; else, detail=info.route; end
fprintf('SMOKE OK: enclosure K=%d M=%d area=%d; boundary=%d; sample=%s (%s)\n', ...
    containsK,containsM,containsArea,boundaryContains,verdict,detail);
end
