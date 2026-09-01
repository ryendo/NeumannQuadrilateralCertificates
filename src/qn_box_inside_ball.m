function tf = qn_box_inside_ball(center,half_widths,radius)
% First discard condition in (46): the complete box lies in the closed ball.

C=qn_global_constants();
if nargin<3, radius=C.rho_sharp_over_2; end
p=qn_interval_box(center,half_widths);
r2=intval('0');
for k=1:4, r2=r2+intval(sup(abs(p(k))))^2; end
tf = sup(sqrt(r2)+C.ulp_pad) <= inf(radius);
end
