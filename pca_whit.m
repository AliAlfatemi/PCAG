function [Out]=pca_whit(A,n)
sum_explained = 0;
idx = 0;
X = zscore(A);
%scatter(X(:,10000),X(:,11000));
[coeff, score, latent, tsquared, explained, mu] = pca(X);
%explained = latent/sum(latent) * 100;
%pareto(latent)
Out=score(:, 1:n)*coeff(:, 1:n)';
%scatter(Out(:,10000),Out(:,11000));



%biplot(V(:,1:121),'scores',U(:,1:121),'varlabels',{'v_1','v_2','v_3','v_4'});
