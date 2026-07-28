function report = qn_global_vectorization_test(repo_root)
% Regression test for the vectorized certified global assembly.

if nargin<1, repo_root=fileparts(fileparts(mfilename('fullpath'))); end
addpath(fullfile(repo_root,'src','common'));
addpath(fullfile(repo_root,'src','global'));

centers={ ...
    zeros(4,1), ...
    [0.10;0.05;0.04;0.02], ...
    [0.12957;0.56535;0.23984;-0.01376], ...
    [-0.17123475299635676;0.20062825348897717;0.0;-0.022264254417977088]};
halfwidths={zeros(4,1),0.01*ones(4,1), ...
    [0.008638;0.017132;0.034263;0.013759],1e-6*ones(4,1)};
% Put the last center exactly on the c_3=0 triangle face.
a=centers{4}(1); b=centers{4}(2); d=centers{4}(4);
centers{4}(3)=(1-a*a-d*d-b*(1-d+a))/(1-d-a);

t0=tic;
for qcase=1:numel(centers)
    center=centers{qcase}; hw=halfwidths{qcase};
    [K,M,area]=qn_km_enclosure(center,hw);
    [Kf,Mf]=qn_km_float(center);
    areaAtCenter=1-center(1)^2-center(4)^2;
    assert(all(all(in(Kf,K))) && all(all(in(Mf,M))) && in(areaAtCenter,area), ...
        'Vectorized certified enclosure misses its center evaluation.');
end
elapsedSeconds=toc(t0);

% The last box genuinely crosses the degenerate triangle face.
p=qn_interval_box(centers{4},halfwidths{4});
c3=intval('1')-p(1)^2-p(4)^2-(p(2)+p(3))*(intval('1')-p(4)) ...
    -p(1)*(p(2)-p(3));
assert(inf(c3)<=0 && sup(c3)>=0,'Boundary regression box does not cross c_3=0.');
[boundaryVerdict,boundaryInfo]=qn_certify_box(centers{4},halfwidths{4});
assert(strcmp(boundaryVerdict,'cert'),'Vectorized boundary box was not certified.');

report=struct('cases',numel(centers),'all_centers_contained',true, ...
    'boundary_verdict',boundaryVerdict, ...
    'boundary_info',boundaryInfo,'elapsed_seconds',elapsedSeconds);
fprintf('GLOBAL VECTORIZATION OK: cases=%d elapsed=%.3fs boundary=%s\n', ...
    report.cases,elapsedSeconds,boundaryVerdict);
end
