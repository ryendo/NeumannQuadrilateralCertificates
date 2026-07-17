function report=dq2_remainder_vectorization_test(repo_root)
% Compare vectorized Taylor-cell bounds with the scalar certified reference.

if nargin<1, repo_root=fileparts(fileparts(mfilename('fullpath'))); end
addpath(fullfile(repo_root,'src','local'));
PI=intval('pi'); rho=intval('3232')/(intval('27')*PI^6);
x=[infsup(-0.01,0.01);infsup(-0.01,0.01);infsup(-0.01,0.01)];
nrm=sqrt(intval('1')+x(1)^2+x(2)^2+x(3)^2);
E=[intval('1')/nrm;x(1)/nrm;x(2)/nrm;x(3)/nrm];

t=tic; [Ko,Mo,bo]=dq2_bound_taylor_remainder(E,rho,2,10); oldSeconds=toc(t);
t=tic; [Kn,Mn,bn]=dq2_bound_taylor_remainder_vectorized(E,rho,2,10); newSeconds=toc(t);
old=[Ko(:);Mo(:);bo(:)]; new=[Kn(:);Mn(:);bn(:)];
assert(all(isfinite(new)) && all(new>=0),'Vectorized remainder radii are invalid.');
assert(all(new<=old),'Vectorized cell evaluation did not refine the scalar reference.');
report=struct('componentwise_subset',true,'old_seconds',oldSeconds, ...
    'new_seconds',newSeconds,'speedup',oldSeconds/newSeconds, ...
    'largest_relative_reduction',max((old-new)./max(1,old)));
fprintf('LOCAL REMAINDER VECTORIZATION OK: %.3fx (%.3fs -> %.3fs), reduction %.3g\n', ...
    report.speedup,oldSeconds,newSeconds,report.largest_relative_reduction);
end
