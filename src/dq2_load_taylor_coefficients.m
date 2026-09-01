function TC = dq2_load_taylor_coefficients()
% Build dense interval coefficient matrices from taylor_coefficients.mat.
% TC.MEXP: 693x4 exponents. TC.CK{j+1}: 15x693 intval, etc. (j = 0..9)
%
% PDF notation:
%   Taylor expansions of K(p), M(p)                (17)
%   Taylor quotients entering C_t(nu)              (27)
%
% Evaluation:
%   taylor_coefficients.mat stores symbolic monomial coefficients generated
%   from the exact trigonometric kernels.  This loader converts every rational
%   coefficient and pi-power into INTLAB intervals.  No floating-point
%   coefficient is trusted as exact after this point.
persistent CACHE
if ~isempty(CACHE), TC = CACHE; return; end
repo_root = fileparts(fileparts(mfilename('fullpath')));
table_path=fullfile(repo_root, 'data', 'taylor_coefficients.mat');
expected_sha256='c753842c1a98e6b7b9c811968755b89bcc043f48645798487cfb2a48b14bb50e';
actual_sha256=qn_sha256_file(table_path);
assert(strcmp(actual_sha256,expected_sha256), ...
    'Taylor coefficient table SHA-256 mismatch.');
S = load(table_path);
PI = intval('pi');
M = size(S.MEXP, 1);
TC.MEXP = S.MEXP;
TC.J = 9;
TC.CK = cell(10,1); TC.CM = cell(10,1); TC.Cb = cell(10,1);
for j = 0:9
    TC.CK{j+1} = intval(zeros(15, M));
    TC.CM{j+1} = intval(zeros(15, M));
    TC.Cb{j+1} = intval(zeros(5, M));
end
% pi-power cache
pmin = min(S.TERMS(:,7)); pmax = max(S.TERMS(:,7));
pip = containers.Map('KeyType','double','ValueType','any');
for p = pmin:2:pmax
    pip(p) = PI^p;
end
for r = 1:size(S.TERMS,1)
    fam = S.TERMS(r,1); ent = S.TERMS(r,2); j = S.TERMS(r,3);
    mi = S.TERMS(r,4); num = S.TERMS(r,5); den = S.TERMS(r,6); pw = S.TERMS(r,7);
    val = (intval(num)/intval(den))*pip(pw);  % symbolic rational*pi^pw -> INTLAB interval
    switch fam
        case 1, TC.CK{j+1}(ent, mi) = TC.CK{j+1}(ent, mi) + val; % K_j(e)
        case 2, TC.CM{j+1}(ent, mi) = TC.CM{j+1}(ent, mi) + val; % M_j(e)
        case 3, TC.Cb{j+1}(ent, mi) = TC.Cb{j+1}(ent, mi) + val; % b_j(e)
    end
end
pr = zeros(15,2); c = 0;
for i = 1:5
    for j2 = i:5
        c = c + 1; pr(c,:) = [i j2];
    end
end
TC.pr = pr;
CACHE = TC;
end
