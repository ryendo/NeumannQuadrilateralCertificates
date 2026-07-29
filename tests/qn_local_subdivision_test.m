function report=qn_local_subdivision_test()
% Regression test for mixed S-only and veigs-undecided radial levels.

outdir=tempname;
mkdir(outdir);
cleanup=onCleanup(@() rmdir(outdir,'s')); %#ok<NASGU>
qn_local_certificate_cover(20,13824,12,outdir);
T=readtable(fullfile(outdir,'res_020.csv'));
assert(height(T)==1 && T.box_id==20 && T.ok==1, ...
    'The representative mixed-level box was not certified.');
assert(T.max_subdivision_depth>=1, ...
    'The regression box did not exercise local subdivision.');
assert(isfile(fullfile(outdir,'done_020.txt')), ...
    'The local worker did not reach its completion marker.');
report=struct('box_id',T.box_id,'infS',T.inf_S, ...
    'subdivision_depth',T.max_subdivision_depth);
fprintf('LOCAL SUBDIVISION OK: box=%d S=%.17g depth=%d\n', ...
    report.box_id,report.infS,report.subdivision_depth);
end
