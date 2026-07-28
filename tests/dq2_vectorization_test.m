function report=dq2_vectorization_test(repo_root)
% Compare vectorized local polynomial evaluation with the scalar reference.

if nargin<1, repo_root=fileparts(fileparts(mfilename('fullpath'))); end
addpath(fullfile(repo_root,'src','local'));
TC=dq2_load_taylor_coefficients();
xI=[infsup(-0.01,0.01);infsup(-0.01,0.01);infsup(-0.01,0.01)];

EI=dq2_face_direction(xI,1,1);
t=tic; oldI=eval_old(TC,EI); oldIntervalSeconds=toc(t);
t=tic; newI=eval_new(TC,EI); newIntervalSeconds=toc(t);
compare_interval_families(oldI,newI);

EP=dq2_face_direction(hessianinit(intval(mid(xI))),1,1);
t=tic; oldP=eval_old(TC,EP); oldPointHessianSeconds=toc(t);
t=tic; newP=eval_new(TC,EP); newPointHessianSeconds=toc(t);
compare_hessian_families(oldP,newP);

EG=dq2_face_direction(hessianinit(xI),1,1);
t=tic; oldG=eval_old(TC,EG); oldBoxHessianSeconds=toc(t);
t=tic; newG=eval_new(TC,EG); newBoxHessianSeconds=toc(t);
compare_hessian_families(oldG,newG);

oldTotal=oldIntervalSeconds+oldPointHessianSeconds+oldBoxHessianSeconds;
newTotal=newIntervalSeconds+newPointHessianSeconds+newBoxHessianSeconds;
report=struct('all_interval_overlaps',true,'old_seconds',oldTotal, ...
    'new_seconds',newTotal,'speedup',oldTotal/newTotal);
fprintf('LOCAL VECTORIZATION OK: %.3fx (%.3fs -> %.3fs)\n', ...
    report.speedup,oldTotal,newTotal);
end

function out=eval_old(TC,E)
[out.K,out.M,out.b]=dq2_evaluate_taylor_coefficients(TC,E(1),E(2),E(3),E(4));
end

function out=eval_new(TC,E)
[out.K,out.M,out.b]=dq2_evaluate_taylor_coefficients_vectorized(TC,E(1),E(2),E(3),E(4));
end

function compare_interval_families(a,b)
names={'K','M','b'};
for f=1:numel(names)
    x=a.(names{f}); y=b.(names{f});
    for j=1:10
        assert(overlap(x{j},y{j}),'Local interval polynomial values do not overlap.');
    end
end
end

function compare_hessian_families(a,b)
names={'K','M','b'};
for f=1:numel(names)
    x=a.(names{f}); y=b.(names{f});
    for j=1:10
        assert(overlap(x{j}.x,y{j}.x),'Local Hessian values do not overlap.');
        assert(overlap(x{j}.dx,y{j}.dx),'Local gradients do not overlap.');
        assert(overlap(x{j}.hx,y{j}.hx),'Local Hessians do not overlap.');
    end
end
end

function tf=overlap(x,y)
xi=inf(x); xs=sup(x); yi=inf(y); ys=sup(y);
tf=all(xi(:)<=ys(:)) && all(yi(:)<=xs(:));
end
