function [left,right] = qn_bisect_box(box,k)
% Source bisection. qn_interval_box outward-rounds both child descriptors.

if nargin<2 || isempty(k), [~,k]=max(box.half_widths); end
old=box.half_widths(k); h=old/2;
left=box; right=box; left.depth=box.depth+1; right.depth=box.depth+1;
left.half_widths(k)=h; right.half_widths(k)=h;
left.center(k)=box.center(k)-h; right.center(k)=box.center(k)+h;
end
