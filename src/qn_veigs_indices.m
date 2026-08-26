function [selected, info] = qn_veigs_indices(A, B, requested_indices)
% Verified enclosures of selected low-end generalized eigenvalue indices.
%
% Calls veigs for the first max(requested_indices) eigenvalues and then uses
% the returned certified index data to select the requested enclosures.  A
% scalar bound representing a complete returned cluster is accepted for each
% requested index in that cluster; this is rigorous, although it may be too
% wide for a later strict inequality.

if exist('veigs','file') ~= 2
    error('qn:VEIGSUnavailable', ...
        'veigs is not on the MATLAB path. Pass its root to QuadrilateralProofRunner.');
end
if ~isequal(size(A),size(B)) || size(A,1) ~= size(A,2)
    error('qn:VEIGSInput','K and M must be square matrices of equal size.');
end

requested_indices = double(requested_indices(:)');
if isempty(requested_indices) || any(~isfinite(requested_indices)) || ...
        any(requested_indices ~= floor(requested_indices)) || ...
        any(requested_indices < 1) || any(requested_indices > size(A,1)) || ...
        numel(unique(requested_indices)) ~= numel(requested_indices)
    error('qn:VEIGSInput','Requested indices must be distinct integers between 1 and matrix size.');
end

% The mathematical pencils are symmetric. Hull with the transpose makes
% that symmetry exact at the interval-data level without discarding either
% enclosure.
A = hull(intval(A),intval(A)');
B = hull(intval(B),intval(B)');

try
    [bounds, indices] = veigs(A,B,max(requested_indices),'sa');
catch ME
    error('qn:VEIGSFailure','veigs failed: %s',ME.message);
end

bounds = bounds(:);
indices = double(indices(:)');
if isempty(indices) || any(~isfinite(indices)) || ...
        any(indices ~= floor(indices)) || any(indices < 1) || ...
        any(indices > size(A,1)) || numel(unique(indices)) ~= numel(indices)
    error('qn:VEIGSIndexData','veigs returned invalid certified index data.');
end
if isempty(bounds) || any(~isfinite(inf(bounds))) || ...
        any(~isfinite(sup(bounds))) || any(inf(bounds) > sup(bounds))
    error('qn:VEIGSBounds','veigs returned invalid eigenvalue bounds.');
end
selected = intval(zeros(1,numel(requested_indices)));
if numel(bounds) == numel(indices)
    for j = 1:numel(requested_indices)
        position = find(indices == requested_indices(j),1);
        if isempty(position)
            error('qn:VEIGSIndex', ...
                'veigs did not certify requested eigenvalue index %d.',requested_indices(j));
        end
        selected(j) = bounds(position);
    end
elseif numel(bounds) == 1
    % One interval may represent the complete returned cluster.
    for j = 1:numel(requested_indices)
        if ~any(indices == requested_indices(j))
            error('qn:VEIGSIndex', ...
                'veigs did not certify requested eigenvalue index %d.',requested_indices(j));
        end
        selected(j) = bounds(1);
    end
else
    error('qn:VEIGSShape', ...
        'veigs returned %d bounds for %d indices.',numel(bounds),numel(indices));
end

info = struct('requested_indices',requested_indices,'indices',indices, ...
    'bounds_inf',inf(bounds(:))','bounds_sup',sup(bounds(:))','commit', ...
    '6556d39a0d9819bb172d232062b698aa76e420f6');
end
