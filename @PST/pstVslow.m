function [V_s, lambda, MK, K, M] = pstVslow(pst, nc)

assert(nc>1);
[MK, K, M] = PST.svm_em(pst);

%Mg = full(diag(M)); invMsqrt = diag(1./sqrt(Mg));
%MKM = diag(sqrt(Mg))*MK*invMsqrt; % K is non-symmetric, so MKM won't be symmetric anyway
%[eig_vec, lambda] = eig(full(MKM));
[eig_vec, lambda]  = eig(full(MK));
lambda = diag(lambda);
[~, k] = sort(real(lambda), 'descend');
lambda = lambda(k);

V_s = eig_vec(:, k(1:nc));  % keep the first nc eigenvectors
V_s = real(V_s);  % just in case
lambda = lambda(1:nc);

if ~all(abs(imag(lambda))<1)
lambda = real(lambda);
return;
end

assert(all(abs(imag(lambda))<1));  %sometimes not all computed eigenvalues are purely real...
lambda = real(lambda);
end

