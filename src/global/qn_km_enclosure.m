function [K,M,area,info] = qn_km_enclosure(center,half_widths)
% Certified first-order mean-value enclosure of the 5x5 pencil.
%
% Df(B) is evaluated by INTLAB gradient arithmetic on the entire box, hence
% f(c)+Df(B)(B-c) contains the complete parameter dependence without a
% Taylor remainder. The spatial GL truncation is padded separately.

center=center(:); half_widths=half_widths(:); pbox=qn_interval_box(center,half_widths);
pc=intval(zeros(4,1));
for k=1:4, pc(k)=intval(sprintf('%.17g',center(k))); end
[Kc,Rc,mc]=qn_assemble_ad(pc);
[Kg,Rg,mg]=qn_assemble_ad(gradientinit(pbox));
[epsK,epsM,padinfo]=qn_gl_pad(pbox,20);

K=intval(zeros(5)); R=intval(zeros(5)); means=intval(zeros(5,1));
padK=infsup(-sup(epsK),sup(epsK)); padM=infsup(-sup(epsM),sup(epsM));
for i=1:5
    means(i)=mean_value(mc{i},mg{i},half_widths)+padM;
    for j=1:5
        K(i,j)=mean_value(Kc{i,j},Kg{i,j},half_widths)+padK;
        R(i,j)=mean_value(Rc{i,j},Rg{i,j},half_widths)+padM;
    end
end
area=intval('1')-pbox(1)^2-pbox(4)^2;
if inf(area)<=0, error('qn:Area','Area is not certifiably positive.'); end
M=R-(means*means')/area;
info=padinfo;
end

function y=mean_value(at_center,on_box,hw)
y=at_center; g=on_box.dx;
for k=1:4
    h=intval(sprintf('%.17g',hw(k)));
    y=y+g(k)*infsup(-sup(h),sup(h));
end
end
