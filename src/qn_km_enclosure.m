function [K,M,q_box,info] = qn_km_enclosure(center,half_widths)
% Certified first-order mean-value enclosure of the 5x5 pencil.
% The output q_box is [q]_B in (44), where q(p)=|Q_p|=1-a^2-d^2.
% Its endpoints are the paper's lower(q_B) and upper(q_B).  The symbol A is
% reserved for the pulled-back form matrix introduced in Section 4.3.
%
% Df(B) is evaluated by INTLAB gradient arithmetic on the entire box, hence
% f(c)+Df(B)(B-c) contains the complete parameter dependence without a
% Taylor remainder. Spatial GL truncation pads are applied both to the
% center values and to the box-uniform gradients used by the mean-value form.

center=center(:); half_widths=half_widths(:); pbox=qn_interval_box(center,half_widths);
pc=intval(zeros(4,1));
for k=1:4, pc(k)=intval(sprintf('%.17g',center(k))); end
[Kc,Rc,mc]=qn_assemble_interval(pc);
[Kg,Rg,mg]=qn_assemble_interval(gradientinit(pbox));
[epsK,epsM,padinfo,epsGradK,epsGradM]=qn_gl_pad(pbox,20);

K=intval(zeros(5)); R=intval(zeros(5)); means=intval(zeros(5,1));
for i=1:5
    means(i)=mean_value(mc{i},mg{i},half_widths,epsM,epsGradM);
    for j=1:5
        K(i,j)=mean_value(Kc{i,j},Kg{i,j},half_widths,epsK,epsGradK);
        R(i,j)=mean_value(Rc{i,j},Rg{i,j},half_widths,epsM,epsGradM);
    end
end
q_box=intval('1')-pbox(1)^2-pbox(4)^2;
if ~isfinite(inf(q_box)) || ~isfinite(sup(q_box)) || inf(q_box)<=0
    error('qn:QNotPositive', ...
        'The condition inf([q]_B)>0 in (44) is not certified.');
end
M=R-(means*means')/q_box;
info=padinfo;
end

function y=mean_value(at_center,on_box,hw,epsValue,epsGradient)
padValue=infsup(-sup(epsValue),sup(epsValue));
padGradient=infsup(-sup(epsGradient),sup(epsGradient));
y=at_center+padValue; g=on_box.dx;
for k=1:4
    h=intval(sprintf('%.17g',hw(k)));
    y=y+(g(k)+padGradient)*infsup(-sup(h),sup(h));
end
end
