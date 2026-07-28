function [lambda1, info] = qn_veigs_smallest(A, B)
% Verified enclosure of the index-1 generalized eigenvalue using veigs.
%
% The dependency is pinned by the repository documentation to
% yuuka-math/veigs commit 6556d39a0d9819bb172d232062b698aa76e420f6.
% veigs may return a cluster containing more indices than requested.  This
% wrapper accepts that case only when the returned data unambiguously
% contains index 1.

if exist('veigs','file') ~= 2
    error('qn:VEIGSUnavailable', ...
        'veigs is not on the MATLAB path. Pass its root to QuadrilateralProofRunner.');
end
if ~isequal(size(A),size(B)) || size(A,1) ~= size(A,2)
    error('qn:VEIGSInput','K and M must be square matrices of equal size.');
end

% The mathematical pencils are symmetric.  Hull with the transpose makes
% that symmetry exact at the interval-data level without discarding either
% enclosure.
A = hull(intval(A),intval(A)');
B = hull(intval(B),intval(B)');

try
    [bounds, indices] = veigs(A,B,1,'sa');
catch ME
    error('qn:VEIGSFailure','veigs failed: %s',ME.message);
end

indices = double(indices(:)');
position = find(indices == 1,1);
if isempty(position)
    error('qn:VEIGSIndex','veigs did not certify a cluster containing index 1.');
end
if numel(bounds) == numel(indices)
    lambda1 = bounds(position);
elseif numel(bounds) == 1
    % A single interval may represent the complete returned cluster.
    lambda1 = bounds(1);
else
    error('qn:VEIGSShape', ...
        'veigs returned %d bounds for %d indices.',numel(bounds),numel(indices));
end
if ~isfinite(inf(lambda1)) || ~isfinite(sup(lambda1)) || inf(lambda1) > sup(lambda1)
    error('qn:VEIGSBounds','veigs returned an invalid enclosure for index 1.');
end

info = struct('indices',indices,'bounds_inf',inf(bounds(:))', ...
    'bounds_sup',sup(bounds(:))','commit', ...
    '6556d39a0d9819bb172d232062b698aa76e420f6');
end
