function report=dq2_remainder_vectorization_test(repo_root)
% Check validity and subdivision monotonicity of Taylor-cell bounds.

if nargin<1, repo_root=fileparts(fileparts(mfilename('fullpath'))); end
addpath(fullfile(repo_root,'src','local'));
PI=intval('pi'); rho=intval('3232')/(intval('27')*PI^6);
x=[infsup(-0.01,0.01);infsup(-0.01,0.01);infsup(-0.01,0.01)];
E=dq2_face_direction(x,1,1);

t=tic; [Kc,Mc,bc]=dq2_bound_taylor_remainder_vectorized(E,rho,1,10); coarseSeconds=toc(t);
t=tic; [Kf,Mf,bf]=dq2_bound_taylor_remainder_vectorized(E,rho,2,10); fineSeconds=toc(t);
coarse=[Kc(:);Mc(:);bc(:)]; fine=[Kf(:);Mf(:);bf(:)];
assert(all(isfinite(fine)) && all(fine>=0),'Taylor remainder radii are invalid.');
tolerance=64*eps(max(1,coarse));
assert(all(fine<=coarse+tolerance), ...
    'Taylor-cell subdivision unexpectedly enlarged a remainder bound.');
report=struct('valid_nonnegative_bounds',true,'subdivision_monotone',true, ...
    'coarse_seconds',coarseSeconds,'fine_seconds',fineSeconds, ...
    'largest_relative_reduction',max((coarse-fine)./max(1,coarse)));
fprintf('LOCAL REMAINDER OK: subdivision monotone (%.3fs -> %.3fs), reduction %.3g\n', ...
    coarseSeconds,fineSeconds,report.largest_relative_reduction);
end
