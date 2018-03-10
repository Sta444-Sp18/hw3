library(magrittr)
library(fields)

set.seed(12413344)

x = seq(0,3, length.out = 100)

d = dist(x) %>% as.matrix()

sigma2=0.5
l = 3.5
sigma2_w = 0.1

beta0=5
beta1=-0.5

cov = sigma2*exp(-(d*l)^2) + ifelse(d==0,1e-6,0)

mu = beta0 + beta1 * x

y = mu + t(chol(cov)) %*% rnorm(length(x))
y2 = y +rnorm(length(x),0,sigma2_w)




x_pred = c(0.01,2.999,rnorm(20,0.5,0.33), rnorm(20, 2.2,0.33), runif(20, 0, 3)) %>% sort() %>% .[. >= 0]
d_pred = dist(x_pred) %>% as.matrix()
cov_pred = sigma2*exp(-(d_pred*l)^2) + ifelse(d_pred==0,sigma2_w,0)

d_o_pred = rdist(x, x_pred) %>% as.matrix()
cov_o_pred = sigma2*exp(-(d_o_pred*l)^2) + ifelse(d_o_pred==0,sigma2_w,0)
cov_pred_o = t(cov_o_pred)

cond_mu  = (beta0 + beta1 * x_pred) + cov_pred_o %*% solve(cov) %*% (y - mu)
cond_cov = cov_pred - cov_pred_o %*% solve(cov) %*% cov_o_pred

y_pred = cond_mu + t(chol(cond_cov)) %*% rnorm(length(x_pred))


plot(x,y,type='l')
points(x_pred, y_pred, pch=16, col='red')

gp = data.frame(x = x_pred, y = y_pred)
gp_truth = data.frame(x = x, y = y)

saveRDS(gp,"gp.rds")
saveRDS(gp_truth,"gp_truth.rds")
