# response variable first, predictor variables following
fishdata = data.frame(read.delim("Data/lipidoxidation.txt", header = TRUE, sep = "",
                                 dec = ".", fill = TRUE, comment.char = "",
                                 stringsAsFactors=FALSE))
# number of predictors
p = ncol(fishdata) - 1
# number of samples
n = nrow(fishdata)
scaled.fish = scale(fishdata)

# Cross-validation to find delta
deltas = seq(0, 0.02, length.out=(p+1))
CV = c()
for (delta in deltas) {
  CV.temp = c()
  for (k in seq(K)) {
    holdout = data.matrix(scaled.fish[folds[[k]],])
    holdoutZ = holdout[,2:(p+1)]
    
    # training data
    tr_ = data.matrix(scaled.fish[-folds[[k]],])
    tr.Z = tr_[,2:(p+1)]
    
    prod_ = t(tr.Z)%*%tr.Z + nrow(tr.Z)*delta*diag(p)
    # parameters of regression model
    beta.rr = solve(prod_, tol=0)%*%t(tr.Z)%*%tr_[,1]
    
    # test and find goodness of training model on holdout
    temp.pred = holdoutZ%*%beta.rr
    CV_ = t(holdout[,1] - temp.pred)%*%(holdout[,1] - temp.pred)
    CV.temp = c(CV.temp, CV_)
  }
  CV = c(CV, mean(CV.temp))
}
best_delta = deltas[which.min(CV)]
print(best_delta)
# find minimal ssr
plot(deltas, CV, ylim=c(0,10))

# regression
rr.beta = solve(t(Z)%*%Z + nrow(Z)*best_delta*diag(p), tol=0)%*%t(Z)%*%y
rr.y = Z%*%rr.beta

# plot
plot.new()
plot.window(xlim = c(min(rr.y)-.1, max(rr.y)+.1),
            ylim = c(min(scaled.fish[,1])-.1, max(scaled.fish[,1]))+.1)
title(main='Plot of the TBARS-values PCR predicted TBARS',
      xlab = 'Predicted TBARS', ylab = 'TBARS')
axis(1)
axis(2)
lines(rr.y, scaled.fish[,1], 'p')
lines(c(-10, 10), c(-10, 10), col='red')
