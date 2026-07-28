function [U,V,W] = qn_gl_tensor20()
% Cached interval tensor grid for the certified 20-point GL assembly.

% qn_gauss_legendre_20 encloses every node and weight.  Reusing the same
% intval arrays changes neither the quadrature rule nor its enclosure.

persistent U_cached V_cached W_cached
if isempty(U_cached)
    [x,w]=qn_gauss_legendre_20(); n=size(x,1);
    % INTLAB's repmat and colon indexing compact thin zero/constant arrays in
    % some releases.  Repeat the directed endpoints and reconstruct intval
    % arrays explicitly so the tensor dimensions are preserved.
    Ug=infsup(repmat(inf(x),1,n),repmat(sup(x),1,n));
    Vg=infsup(repmat(inf(x).',n,1),repmat(sup(x).',n,1));
    Wu=infsup(repmat(inf(w),1,n),repmat(sup(w),1,n));
    Wv=infsup(repmat(inf(w).',n,1),repmat(sup(w).',n,1));
    Wg=Wu.*Wv;
    U_cached=infsup(reshape(inf(Ug),[],1),reshape(sup(Ug),[],1));
    V_cached=infsup(reshape(inf(Vg),[],1),reshape(sup(Vg),[],1));
    W_cached=infsup(reshape(inf(Wg),[],1),reshape(sup(Wg),[],1));
end
U=U_cached; V=V_cached; W=W_cached;
end
