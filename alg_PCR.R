# This is a translation of https://towardsdatascience.com/implementing-pca-from-scratch-fb434f1acbaa
# into R code. Formerly, this was python code. A reference class is used.

library(dplyr)
library(Matrix)
setClassUnion("df_type", c("matrix", "data.frame", "dgCMatrix"))
PCA = setRefClass("PCA", fields = list(data="df_type"),
            methods=list(
              
              standardize_data = function(X=data){
                apply(X, 2, function(x) (x - mean(x)) / sd(x))
              },
              
              get_covariance_matrix = function(X=data){
                cov(X)
              },
              
              get_eigenvectors = function(X=data, n_cols){
                eigenvalues = eigen(X)$values
                eigenvectors = eigen(X)$vectors
                cols = sort(eigenvalues, decreasing=TRUE, index.return=TRUE)$ix
                eigenvectors[,cols[1:n_cols]]
              },
              
              project_matrix = function(X=data, eigenvectors){
                X %*% eigenvectors
              }
            ))

fit_transform = function(X, n_cols){
  eigenvectors = X %>% (apply_PCA$standardize_data) %>% (apply_PCA$get_covariance_matrix) %>% (function(x) apply_PCA$get_eigenvectors(x, n_cols))
  apply_PCA$project_matrix(as.matrix(X), eigenvectors) 
}

data(iris)
Y = iris[,ncol(iris)]
X = iris[,-ncol(iris)]
plot(X[,c(1,2)],col=Y)

apply_PCA = PCA(data=as.matrix(iris))
plot(fit_transform(X, n_cols=2),col=Y)
