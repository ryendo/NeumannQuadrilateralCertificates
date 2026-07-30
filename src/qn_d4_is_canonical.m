function tf = qn_d4_is_canonical(center)
% Source-code convention for one representative of each D4 orbit.

c=center(:);
a=c(1); b=c(2); cc=c(3); d=c(4);
images=[a, a,  a,  a, -a, -a, -a, -a; ...
        b, b, -b, -b, -cc,-cc, cc, cc; ...
        cc,-cc,cc,-cc, -b,  b, -b,  b; ...
        d,-d, -d,  d,   d, -d, -d,  d];
keys=round(images.'*1e12);
best=keys(1,:);
for k=2:size(keys,1)
    first_difference=find(keys(k,:)~=best,1);
    if ~isempty(first_difference) && keys(k,first_difference)>best(first_difference)
        best=keys(k,:);
    end
end
tf=isequal(round(c.'*1e12),best);
end
