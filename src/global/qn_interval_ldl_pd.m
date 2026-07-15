function ok = qn_interval_ldl_pd(A)
% Interval LDL' test without pivoting. Success certifies A(p) positive definite.

n=size(A,1); L=intval(eye(n)); D=intval(zeros(n,1)); ok=true;
for j=1:n
    s=intval('0');
    for k=1:j-1, s=s+L(j,k)^2*D(k); end
    D(j)=A(j,j)-s;
    if inf(D(j))<=0, ok=false; return; end
    for i=j+1:n
        s=intval('0');
        for k=1:j-1, s=s+L(i,k)*L(j,k)*D(k); end
        L(i,j)=(A(i,j)-s)/D(j);
    end
end
end
