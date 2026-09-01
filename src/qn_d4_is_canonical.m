function tf = qn_d4_is_canonical(center)
% True for one distinct representative of the D_4 action preceding (13).

c=center(:);
images=qn_d4_images(c);
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
