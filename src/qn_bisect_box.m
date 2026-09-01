function [left,right] = qn_bisect_box(box,k)
% Covering bisection corresponding to (47).
% Child descriptors are constructed from the certified parent endpoints so
% that their INTLAB interpretations cover the parent despite binary rounding.

if nargin<2 || isempty(k), [~,k]=max(box.half_widths); end
parent_interval=qn_interval_box(box.center,box.half_widths);
lo=inf(parent_interval(k)); hi=sup(parent_interval(k));
junction=box.center(k);

left=box; right=box;
left.depth=box.depth+1; right.depth=box.depth+1;
[left.center(k),left.half_widths(k)]=covering_descriptor(lo,junction);
[right.center(k),right.half_widths(k)]=covering_descriptor(junction,hi);
end

function [center,half_width]=covering_descriptor(lo,hi)
center=(lo+hi)/2;
half_width=max(center-lo,hi-center);
for attempt=1:8
    center_interval=intval(sprintf('%.17g',center));
    half_width_interval=intval(sprintf('%.17g',half_width));
    box_interval=center_interval+ ...
        infsup(-sup(half_width_interval),sup(half_width_interval));
    if inf(box_interval)<=lo && sup(box_interval)>=hi
        return
    end
    scale=max(abs([lo,hi,center,half_width,1]));
    half_width=half_width+8*eps(scale);
end
error('qn:BisectionRounding', ...
    'Could not construct an outward-rounded child interval for (47).');
end
