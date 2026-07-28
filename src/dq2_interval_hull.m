function X=dq2_interval_hull(a,b)
% Outward-rounded hull of two scalar endpoints.
a=intval(a); b=intval(b);
X=infsup(inf(a),sup(b));
end
