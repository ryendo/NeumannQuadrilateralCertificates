function report = qn_global_vectorization_test(repo_root)
% Regression test for the vectorized certified global assembly.
%
% The legacy scalar jet is retained as an independent interval reference.
% At both point and box inputs, every value and first-derivative enclosure
% produced by the two algebraically identical formulas must overlap.

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

tOldCenter=0; tNewCenter=0; tOldGrad=0; tNewGrad=0;
for qcase=1:numel(centers)
    center=centers{qcase}; hw=halfwidths{qcase};
    pc=point_interval(center); pbox=qn_interval_box(center,hw);

    t=tic; [Kco,Rco,mco]=qn_assemble_ad(pc); tOldCenter=tOldCenter+toc(t);
    t=tic; [Kcn,Rcn,mcn]=qn_assemble_interval_center(pc); tNewCenter=tNewCenter+toc(t);
    assert_center_overlap(Kco,Rco,mco,Kcn,Rcn,mcn);

    t=tic; [Kgo,Rgo,mgo]=qn_assemble_jet1(pbox); tOldGrad=tOldGrad+toc(t);
    t=tic; [Kgn,Rgn,mgn]=qn_assemble_interval_grad(pbox); tNewGrad=tNewGrad+toc(t);
    assert_jet_overlap(Kgo,Rgo,mgo,Kgn,Rgn,mgn);

    [K,M,area]=qn_km_enclosure(center,hw);
    [Kf,Mf]=qn_km_float(center);
    areaAtCenter=intval('1')-pc(1)^2-pc(4)^2;
    assert(all(all(in(Kf,K))) && all(all(in(Mf,M))) && in(mid(areaAtCenter),area), ...
        'Vectorized certified enclosure misses its center evaluation.');
end

% The last box genuinely crosses the degenerate triangle face.
p=qn_interval_box(centers{4},halfwidths{4});
c3=intval('1')-p(1)^2-p(4)^2-(p(2)+p(3))*(intval('1')-p(4)) ...
    -p(1)*(p(2)-p(3));
assert(inf(c3)<=0 && sup(c3)>=0,'Boundary regression box does not cross c_3=0.');
[boundaryVerdict,boundaryInfo]=qn_certify_box(centers{4},halfwidths{4});
assert(strcmp(boundaryVerdict,'cert'),'Vectorized boundary box was not certified.');

report=struct('cases',numel(centers),'all_interval_overlaps',true, ...
    'all_centers_contained',true,'boundary_verdict',boundaryVerdict, ...
    'boundary_info',boundaryInfo,'old_center_seconds',tOldCenter, ...
    'new_center_seconds',tNewCenter,'old_gradient_seconds',tOldGrad, ...
    'new_gradient_seconds',tNewGrad,'center_speedup',tOldCenter/tNewCenter, ...
    'gradient_speedup',tOldGrad/tNewGrad);
fprintf(['GLOBAL VECTORIZATION OK: cases=%d center %.3fx (%.3fs -> %.3fs), ' ...
    'gradient %.3fx (%.3fs -> %.3fs), boundary=%s\n'],report.cases, ...
    report.center_speedup,tOldCenter,tNewCenter,report.gradient_speedup, ...
    tOldGrad,tNewGrad,boundaryVerdict);
end

function p=point_interval(center)
p=intval(zeros(4,1));
for k=1:4, p(k)=intval(sprintf('%.17g',center(k))); end
end

function assert_center_overlap(K1,R1,m1,K2,R2,m2)
for i=1:5
    assert(overlaps(m1{i},m2{i}),'Center mean intervals do not overlap.');
    for j=1:5
        assert(overlaps(K1{i,j},K2{i,j}),'Center stiffness intervals do not overlap.');
        assert(overlaps(R1{i,j},R2{i,j}),'Center raw-mass intervals do not overlap.');
    end
end
end

function assert_jet_overlap(K1,R1,m1,K2,R2,m2)
for i=1:5
    assert(overlaps(m1{i}.v,m2{i}.v),'Mean values do not overlap.');
    for k=1:4
        assert(overlaps(m1{i}.g(k),m2{i}.g(k)),'Mean gradients do not overlap.');
    end
    for j=1:5
        assert(overlaps(K1{i,j}.v,K2{i,j}.v),'Stiffness values do not overlap.');
        assert(overlaps(R1{i,j}.v,R2{i,j}.v),'Raw-mass values do not overlap.');
        for k=1:4
            assert(overlaps(K1{i,j}.g(k),K2{i,j}.g(k)), ...
                'Stiffness gradients do not overlap.');
            assert(overlaps(R1{i,j}.g(k),R2{i,j}.g(k)), ...
                'Raw-mass gradients do not overlap.');
        end
    end
end
end

function tf=overlaps(x,y)
tf=inf(x)<=sup(y) && inf(y)<=sup(x);
end
