function [FR,K0,M0]=dq2_fixed_reference()
% Exact square pencil and symmetry-adapted bases used by Algorithm 1.
PI2=sqr(intval('pi'));
SQ2=sqrt(intval(2));
V5=intval(zeros(5,2)); V5(1,1)=SQ2; V5(2,2)=SQ2;
W5=intval(zeros(5,3)); W5(3,1)=intval(2); W5(4,2)=SQ2; W5(5,3)=SQ2;
D0=intval(diag([PI2,3*PI2,3*PI2]));
FR=struct('V5',V5,'W5',W5,'D0',D0,'I2',intval(eye(2)), ...
    'I3',intval(eye(3)),'PI2',PI2);
if nargout>1, K0=intval(diag([PI2/2,PI2/2,PI2/2,2*PI2,2*PI2])); end
if nargout>2, M0=intval(diag([1/2,1/2,1/4,1/2,1/2])); end
end
