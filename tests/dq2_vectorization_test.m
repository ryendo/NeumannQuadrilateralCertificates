function report=dq2_vectorization_test(repo_root)
% Check that box evaluations contain their center evaluations.

if nargin<1, repo_root=fileparts(fileparts(mfilename('fullpath'))); end
addpath(fullfile(repo_root,'src','local'));
TC=dq2_load_taylor_coefficients();
xI=[infsup(-0.01,0.01);infsup(-0.01,0.01);infsup(-0.01,0.01)];

EI=dq2_face_direction(xI,1,1);
E0=dq2_face_direction(intval(mid(xI)),1,1);
t=tic; boxValues=evaluate(TC,EI); intervalSeconds=toc(t);
centerValues=evaluate(TC,E0);
assert_contains(centerValues,boxValues);

EG=dq2_face_direction(hessianinit(xI),1,1);
EP=dq2_face_direction(hessianinit(intval(mid(xI))),1,1);
t=tic; boxHessian=evaluate(TC,EG); hessianSeconds=toc(t);
centerHessian=evaluate(TC,EP);
assert_hessian_contains(centerHessian,boxHessian);

report=struct('interval_contains_center',true, ...
    'box_hessian_contains_center',true,'interval_seconds',intervalSeconds, ...
    'hessian_seconds',hessianSeconds);
fprintf('LOCAL VECTORIZATION OK: center contained (interval %.3fs, Hessian %.3fs)\n', ...
    intervalSeconds,hessianSeconds);
end

function out=evaluate(TC,E)
[out.K,out.M,out.b]=dq2_evaluate_taylor_coefficients_vectorized(TC,E(1),E(2),E(3),E(4));
end

function assert_contains(point,box)
names={'K','M','b'};
for f=1:numel(names)
    p=point.(names{f}); b=box.(names{f});
    for j=1:numel(p)
        assert(all(in(p{j}(:),b{j}(:))), ...
            'Local interval polynomial evaluation misses its center value.');
    end
end
end

function assert_hessian_contains(point,box)
names={'K','M','b'};
for f=1:numel(names)
    p=point.(names{f}); b=box.(names{f});
    for j=1:numel(p)
        assert(all(in(p{j}.x(:),b{j}.x(:))), ...
            'Local Hessian value enclosure misses its center.');
        assert(all(in(p{j}.dx(:),b{j}.dx(:))), ...
            'Local gradient enclosure misses its center.');
        assert(all(in(p{j}.hx(:),b{j}.hx(:))), ...
            'Local Hessian enclosure misses its center.');
    end
end
end
