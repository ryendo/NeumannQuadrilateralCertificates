function p = qn_interval_box(center, half_widths)
% Convert a floating-point box descriptor into an outward-rounded INTLAB box.

p = intval(zeros(4,1));
for k = 1:4
    c = intval(sprintf('%.17g', center(k)));
    h = intval(sprintf('%.17g', half_widths(k)));
    p(k) = c + infsup(-sup(h), sup(h));
end
end
