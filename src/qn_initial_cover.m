function boxes = qn_initial_cover(n)
% Uniform outward-padded cover of P_C.

C=qn_global_constants(); h=sup(C.pc_half(:)); boxes=cell(n^4,1); r=0;
edges=cell(4,1);
for k=1:4, edges{k}=linspace(-h(k),h(k),n+1); end
for i1=1:n
    for i2=1:n
        for i3=1:n
            for i4=1:n
                r=r+1; idx=[i1 i2 i3 i4]; center=zeros(4,1); hw=zeros(4,1);
                for k=1:4
                    lo=edges{k}(idx(k)); hi=edges{k}(idx(k)+1);
                    center(k)=(lo+hi)/2;
                    hw(k)=(hi-lo)/2 + 8*eps(max(abs([lo hi center(k) 1])));
                end
                boxes{r}=struct('center',center,'half_widths',hw, ...
                    'depth',0,'initial_id',r);
            end
        end
    end
end
end
