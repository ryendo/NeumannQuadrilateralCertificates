function [K,Mraw,means]=qn_assemble_interval(p)
% Vectorized GL assembly for intval or INTLAB gradient parameters.

[U,V,W]=qn_gl_tensor20();
out=qn_quadrilateral_kernels(U,V,p(1),p(2),p(3),p(4));
K=cell(5); Mraw=cell(5); means=cell(5,1); k=0;
for i=1:5
    for j=i:5
        k=k+1;
        K{i,j}=sum(out{k}.*W);
        Mraw{i,j}=sum(out{15+k}.*W);
    end
    means{i}=sum(out{30+i}.*W);
end
for i=1:5
    for j=i+1:5
        K{j,i}=K{i,j}; Mraw{j,i}=Mraw{i,j};
    end
end
end
