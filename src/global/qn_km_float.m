function [K, M, area] = qn_km_float(p)
% Non-certified center assembly. Used only to choose a test frame and split axis.

[xi, wi] = qn_gauss_legendre_20();
x = mid(xi); w = mid(wi);
modes = [1 0; 0 1; 1 1; 2 0; 0 2];
K = zeros(5); Mraw = zeros(5); means = zeros(5,1);
area = 0;
a = p(1); b = p(2); c = p(3); d = p(4);
for iu = 1:20
    u = x(iu);
    for iv = 1:20
        v = x(iv); W = w(iu)*w(iv);
        X = u - a*u - d*v - 2*b*u*v;
        Y = v - d*u + a*v + 2*c*u*v;
        Xu = 1-a-2*b*v; Xv = -(d+2*b*u);
        Yu = -(d-2*c*v); Yv = 1+a+2*c*u;
        J = Xu*Yv-Xv*Yu;
        gvv = Xv^2+Yv^2; guv = Xu*Xv+Yu*Yv; guu = Xu^2+Yu^2;
        area = area + J*W;
        phi = zeros(5,1); phiu = zeros(5,1); phiv = zeros(5,1);
        for q = 1:5
            mm = modes(q,1); nn = modes(q,2);
            cx = cos(mm*pi*(X+0.5)); sx = sin(mm*pi*(X+0.5));
            cy = cos(nn*pi*(Y+0.5)); sy = sin(nn*pi*(Y+0.5));
            phi(q) = cx*cy;
            psiX = -mm*pi*sx*cy; psiY = -nn*pi*cx*sy;
            phiu(q) = psiX*Xu + psiY*Yu;
            phiv(q) = psiX*Xv + psiY*Yv;
        end
        means = means + phi*J*W;
        for i = 1:5
            for j = i:5
                K(i,j) = K(i,j) + (gvv*phiu(i)*phiu(j) ...
                    - guv*(phiu(i)*phiv(j)+phiv(i)*phiu(j)) ...
                    + guu*phiv(i)*phiv(j))/J*W;
                Mraw(i,j) = Mraw(i,j) + phi(i)*phi(j)*J*W;
            end
        end
    end
end
K = K + triu(K,1)'; Mraw = Mraw + triu(Mraw,1)';
M = Mraw - (means*means')/area;
end
