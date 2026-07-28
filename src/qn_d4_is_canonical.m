function tf = qn_d4_is_canonical(center)
% Source-code convention for one representative of each D4 orbit.

c=center(:);
if ~(c(3)>=-1e-12 && c(4)>=-1e-12), tf=false; return; end
img=[-c(1);c(3);c(2);c(4)];
lhs=round([c(1) c(2)]*1e12); rhs=round([img(1) img(2)]*1e12);
tf=(lhs(1)>rhs(1)) || (lhs(1)==rhs(1) && lhs(2)>=rhs(2));
end
