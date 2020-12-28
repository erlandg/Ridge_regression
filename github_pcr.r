# response variable first, predictor variables following
fishdata = data.frame(read.delim("Data/lipidoxidation.txt", header = TRUE, sep = "",
                                 dec = ".", fill = TRUE, comment.char = "",
                                 stringsAsFactors=FALSE))
# number of predictors
p = ncol(fishdata) - 1
# number of samples
n = nrow(fishdata)
scaled.fish = scale(fishdata)

eig.fish = eigen(cov(scaled.fish[,2:100]))
# imported the data, scaled it, and derived its eigenvectors.

Z = data.matrix(scaled.fish[,2:100])
y = scaled.fish[,1]

# K-fold cross validataion
K = 5
folds = split(sample(n, n, replace = FALSE), as.factor(1:K))
K.folds = lapply(folds, function(x) scaled.fish[x, ])

# goodness measure SSR
SSR = c()
for (m in seq(p)) {
  SSR.temp = c()
  for (k in seq(K)) {
    holdout = K.folds[[k]]
    # training data
    tr_ = data.matrix(scaled.fish[-folds[[k]],])
    W.tr = tr_[,2:p+1]%*%eig.fish$vectors[,1:m]
    g.tr = solve(t(W.tr)%*%W.tr, tol=0)%*%t(W.tr)%*%tr_[,1]
    # parameters of training model
    beta.m = eig.fish$vectors[,1:m]%*%g.tr
    # use on holdout
    y.pcr = holdout[,2:p+1]%*%beta.m
    # error measure
    SSR_ = t(holdout[,1] - y.pcr)%*%(holdout[,1] - y.pcr)
    
    SSR.temp = c(SSR.temp, SSR_)
  }
  mean(SSR.temp)
  SSR = c(SSR, mean(SSR.temp))
}

# optimal number of pcs
m = which.min(SSR)
print(m)


# now, for the regression part
W = Z%*%eig.fish$vectors[,1:m]
gamma = solve(t(W)%*%W)%*%t(W)%*%scaled.fish[,1]
pred.y = W%*%gamma


# plot results
plot.new()
plot.window(xlim = c(min(pred.y)-.1, max(pred.y)+.1),
            ylim = c(min(scaled.fish[,1])-.1, max(scaled.fish[,1]))+.1)
title(main='Plot of the TBARS-values PCR predicted TBARS',
      xlab = 'Predicted TBARS', ylab = 'TBARS')
axis(1)
axis(2)
lines(pred.y, scaled.fish[,1], 'p')
lines(c(-10, 10), c(-10, 10), col='red')
