function E=dq2_face_direction(x,axisdim,sgn)
% Map one standard face chart to a direction on S^3.
other=setdiff(1:4,axisdim);
nrm=sqrt(1+x(1)^2+x(2)^2+x(3)^2);
E=cell(4,1);
E{axisdim}=sgn/nrm;
for d=1:3
    E{other(d)}=x(d)/nrm;
end
E=vertcat(E{:});
end
