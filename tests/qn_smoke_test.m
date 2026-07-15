function report = qn_smoke_test(repo_root)
% Quick integration check; requires MATLAB and INTLAB but no full recomputation.

if nargin<1, repo_root=fileparts(fileparts(mfilename('fullpath'))); end
addpath(fullfile(repo_root,'src','local')); addpath(fullfile(repo_root,'src','global'));
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
local=qn_summarize_local_results(fullfile(repo_root,'results','local'));
assert(local.verified,'Saved local certificate is incomplete.');
report=struct('constants',struct('rho_local',mid(C.rho_local),'rho_seam',mid(C.rho_seam)), ...
    'enclosure',struct('K',containsK,'M',containsM,'area',containsArea), ...
    'sample_box',struct('verdict',verdict,'info',info,'gap',gap,'split_dim',sdim), ...
    'local_saved_results_verified',local.verified);
if isfield(info,'reason'), detail=info.reason; else, detail=info.route; end
fprintf('SMOKE OK: enclosure K=%d M=%d area=%d; sample=%s (%s)\n', ...
    containsK,containsM,containsArea,verdict,detail);
end
