function [selected,info]=qn_veigs_indices(A,B,requested)
% Verified bounds for selected generalized-eigenvalue indices.
%
% Each index is targeted separately. Asking veigs for the first k values can
% return one interval for a cluster of k indices, too coarse for an
% index-specific strict inequality.

if exist('veigs','file')~=2
    error('qn:VEIGSUnavailable','veigs is not on the MATLAB path.');
end
if ~isequal(size(A),size(B)) || size(A,1)~=size(A,2)
    error('qn:VEIGSInput','K and M must be square matrices of equal size.');
end
requested=double(requested(:)'); n=size(A,1);
if isempty(requested) || any(requested~=floor(requested)) || ...
        any(requested<1) || any(requested>n) || numel(unique(requested))~=numel(requested)
    error('qn:VEIGSInput','Requested indices must be distinct integers in 1:n.');
end
A=hull(intval(A),intval(A)'); B=hull(intval(B),intval(B)');
approx=sort(real(eig(mid(A),mid(B))));
selected=intval(zeros(1,numel(requested))); certified=[];

for j=1:numel(requested)
    index=requested(j);
    if index==1
        target='sa';
    elseif index==n
        target='la';
    else
        target=approx(index);
    end
    try
        [bounds,indices]=veigs(A,B,1,target);
    catch primary
        try
            [bounds,indices]=veig(A,B,index);
        catch fallback
            error('qn:VEIGSFailure', ...
                'veigs failed for index %d: %s; veig fallback failed: %s', ...
                index,primary.message,fallback.message);
        end
    end
    bounds=bounds(:); indices=double(indices(:)');
    if isempty(indices) || any(indices~=floor(indices)) || ...
            any(indices<1) || any(indices>n) || ~any(indices==index)
        error('qn:VEIGSIndex','veigs did not certify index %d.',index);
    end
    if numel(bounds)==1
        selected(j)=bounds(1);
    elseif numel(bounds)==numel(indices)
        selected(j)=bounds(find(indices==index,1));
    else
        error('qn:VEIGSShape','veigs returned %d bounds for %d indices.', ...
            numel(bounds),numel(indices));
    end
    if ~isfinite(inf(selected(j))) || ~isfinite(sup(selected(j)))
        error('qn:VEIGSBounds','veigs returned a non-finite bound for index %d.',index);
    end
    certified=union(certified,indices);
end

info=struct('requested_indices',requested,'indices',certified,'commit', ...
    '6556d39a0d9819bb172d232062b698aa76e420f6');
end
