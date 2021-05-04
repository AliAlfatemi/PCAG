function [Out]=pca_whit(A,n)
sum_explained = 0;
idx = 0;
X = zscore(A);
%scatter(X(:,10000),X(:,11000));

% De-mean (MATLAB will de-mean inside of PCA, but I want the de-meaned values later)
X = X - mean(X); % Use X = bsxfun(@minus,X,mean(X)) if you have an older version of MATLAB
% Do the PCA
[coeff, score, latent, tsquared, explained, mu] = pca(X);
% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix = cov(X);
[V,D] = eig(covarianceMatrix);
% "coeff" are the principal component vectors.
% These are the eigenvectors of the covariance matrix.
% Compare the columns of coeff and V.
% (Note that the columns are not necessarily in the same *order*,
%  and they might be *lightly different from each other
%  due to floating-point error.)
coeff
V
% Multiply the original data by the principal component vectors
% to get the projections of the original data on the
% principal component vector space. This is also the output "score".
% Compare ...
dataInPrincipalComponentSpace = X*coeff
score
% The columns of X*coeff are orthogonal to each other. This is shown with ...
corrcoef(dataInPrincipalComponentSpace)
% The variances of these vectors are the eigenvalues of the covariance matrix, and are also the output "latent". Compare
% these three outputs
var(dataInPrincipalComponentSpace)'
latent
sort(diag(D),'descend')
Out=score(:, 1:n)*coeff(:, 1:n)';
%scatter(Out(:,10000),Out(:,11000));

%biplot(V(:,1:121),'scores',U(:,1:121),'varlabels',{'v_1','v_2','v_3','v_4'});
