function W = LRSDL_mod_updateW(X, Y_range, Q)

XAXt = (2*X + buildM_2Mbar(X, Y_range, 1))*X';
[V,D] = eig(XAXt,'vector');
[~,I] = sort(D,1,'ascend');
W = V(:,I(1:Q));