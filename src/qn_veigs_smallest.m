function [lambda1, info] = qn_veigs_smallest(A, B)
% Verified enclosure of the index-1 generalized eigenvalue using veigs.
%
% Kept as the global-certificate compatibility wrapper. The local
% certificate calls qn_veigs_indices directly for indices 1 and 3.
[selected,info] = qn_veigs_indices(A,B,1);
lambda1 = selected(1);
end
