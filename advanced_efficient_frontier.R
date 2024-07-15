#install.packages("tseries")
#library(tseries)


rm(list = ls())



#setting the dates
start_date <- as.Date("2019-01-01")
end_date <- as.Date("2020-01-01")



#setting the template of the link
utpl <- "https://query1.finance.yahoo.com/v7/finance/download/%s?period1=%s&period2=%s&interval=1d&events=history&includeAdjustedClose=true"

#setting the period numbers to be inserted in the link
period1 <- as.integer(as.POSIXct(start_date)) 
period2 <-   as.integer(as.POSIXct(end_date))


tickers <- (c("MSFT MCD JPM JNJ") |> strsplit(split= " +"))[[1]] 


filenames <- sapply(tickers,function(x) paste0(x, ".csv"))

data <- "my_data"
if (!dir.exists(data)) {cat("path to the folder created above: ", paste0(getwd(),"/",data)); dir.create(data)} else {"already exists"}


sapply(tickers, function(ticker) download.file(
  sprintf(utpl,ticker,period1,period2), 
  paste0(data, "/", ticker,".csv")))

dir(data)


stocks <- sapply(tickers, function(ticker) as.numeric(read.csv(paste0(data, "/", ticker,".csv"))$'Adj.Close'))


#for the dates we use the microsoft trading days as a reference
rownames(stocks) <- read.csv(paste0(data, "/", "MSFT",".csv"))$'Date'


#downloading CNY/USD
download.file(sprintf(utpl,"CNY=X",period1,period2), paste0(data, "/", "CNY=X",".csv"))
currencies <- (read.csv(paste0(data, "/", "CNY=X",".csv"))$'Adj.Close')
stock_curr_matrix <- cbind(stocks, currencies[-((length(currencies) - nrow(stocks)) : length(currencies))])
colnames(stock_curr_matrix) <- c(tickers, "USD.CNY")



#computing log returns of CNY against USD
stock_curr_matrix <- as.data.frame(stock_curr_matrix)
curr_rets <- diff(log(as.numeric(stock_curr_matrix$USD.CNY)))






#setting the weights
weights <- c(0.25, 0.25, 0.25, 0.25)
benchmark_weights <- c(0.05733366, 0.00583947, 0.01165359, 0.01206661)


rets <- diff(log(stocks))
r <- matrix(colMeans(rets)) # to compute the mean return for each stock
S <- cov(rets)
one <- matrix(rep(1,length(tickers)))




#computing portfolio expected return
ptf_exp_ret <- t(as.matrix(weights)) %*% r

#computing portfolio variance
ptf_var <- t(weights)%*% S %*% weights



#computing the tracking error
tracking_error <- sqrt(t(weights - benchmark_weights)%*% S %*%(weights - benchmark_weights))


#setting parameters to compute optimal weights and minimum variance
s1 <- solve(S)
a <- t(r) %*% s1 %*% r 
b <- t(r) %*% s1 %*% one
c <- t(one) %*% s1 %*% one
mu <- mean(r) # you can choose the value that you prefer
numerator <- t((mu*c - b)%*%t(r) + (a - mu*b)%*%t(one))
denominator <- a * c - b^2


#computing optimal weights (Remember that using this technique you could get negative values as well)
w_star <- s1 %*% (numerator) / as.numeric(denominator)*one

#computing minimum variance
min_var <- t(w_star)%*%S%*%w_star








########Computing  POSITIVE optimal weights #############à
#install.packages("hydroPSO")
library(hydroPSO)


v <- function(w) t(w) %*% S %*% w 
constr1 <- function(w) t(w) %*% one - 1 
constr2 <- function(w) t(w) %*% r- mu 
v.c <- function(w) v(w) + (constr1(w) + constr2(w))^2 

n <- length(tickers)
opt <- hydroPSO(fn = v.c, method="fips", control=list(MinMax="min", write2disk=FALSE), lower= rep(0, n), upper=rep(1, n))
pos_w_star <- opt$par
w_dict <- setNames(pos_w_star, tickers)
sum(pos_w_star)
############################









###### optimal return and variance using positive weights #######
opt_ptf_exp_ret <- t(as.matrix(pos_w_star)) %*% r
opt_ptf_exp_ret

opt_ptf_var <- t(pos_w_star)%*% S %*% pos_w_star
opt_ptf_var
###############################



###generating EFFICIENT FRONTIER ############
n_portfolios <- 10000
w_matrix <- matrix(0, nrow = n_portfolios, ncol = length(tickers))


for (i in 1:n_portfolios) {
  row_sum <- 0
  for (j in 1:(length(tickers) - 1)) {
    w_matrix[i, j] <- runif(1, 0, 1 - row_sum)
    row_sum <- row_sum + w_matrix[i, j]
  }
  w_matrix[i, length(tickers)] <- 1 - row_sum
}

w_matrix <- rbind(w_matrix, pos_w_star)
print(w_matrix)


row_list <- matrix(split(w_matrix, row(w_matrix)), length(w_matrix), 1)
ret.var.matrix <- as.data.frame(t(sapply(row_list, function(w) {c(Returns=t(as.matrix(w)) %*% r, Variances=t(w)%*% S %*% w)})))
plot(sqrt(ret.var.matrix$Variances), ret.var.matrix$Returns, col="orange", xlim = c(0.00065, 0.014), ylim = c(0, 0.0020), xlab = "stdev", ylab = "exp_ret")
points(sqrt(opt_ptf_var), opt_ptf_exp_ret, col="blue", pch=16)
text(sqrt(opt_ptf_var), opt_ptf_exp_ret, pos=3, "P")


#generating a dataframe that indicates var-ret combination for each set of weights
opt.ret.var.matrix <- cbind(row_list, ret.var.matrix)








#plotting just the efficient part
my_data <- ret.var.matrix[order(sqrt(ret.var.matrix$Variances)),]
points(sqrt(my_data[which(my_data$Returns == cummax(my_data$Returns)), ]$Variances), my_data[which(my_data$Returns == cummax(my_data$Returns)), ]$Returns, col="green")
lines(sqrt(my_data[which(my_data$Returns == cummax(my_data$Returns)), ]$Variances), my_data[which(my_data$Returns == cummax(my_data$Returns)), ]$Returns, col="black")

#plotting inefficient part
points(sqrt(ret.var.matrix[ret.var.matrix$Returns <= as.numeric(opt_ptf_exp_ret), ]$Variances), ret.var.matrix[(ret.var.matrix$Returns <= as.numeric(opt_ptf_exp_ret)), ]$Returns, col="grey")







##plotting CAL for minimum-variance-portfolio ####
rf <- .00065
sharpe <- (opt_ptf_exp_ret - rf) / sqrt(opt_ptf_var)   #for each small increment of variance, excess return increases by sharpe_ratio  
abline(rf, sharpe)  
###################




####plotting and computing all the sharpe ratios#####################
sharpes_list <- vector("numeric", nrow(ret.var.matrix))

#sapply(1:nrow(ret.var.matrix), function(i) {
 #abline(rf, (ret.var.matrix[i, ]$Returns - rf) / sqrt(ret.var.matrix[i, ]$Variances))})

for (i in 1:nrow(ret.var.matrix)) {
  sharpes_list[i] <- (ret.var.matrix[i, ]$Returns - rf) / sqrt(ret.var.matrix[i, ]$Variances)
} 

opt_risky_ptf <- cbind(sharpes_list, ret.var.matrix)[which.max(cbind(sharpes_list, ret.var.matrix)$sharpes_list), ]
points(sqrt(opt_risky_ptf$Variances), opt_risky_ptf$Returns, col="red", cex=1, pch=16)
text(sqrt(opt_risky_ptf$Variances), opt_risky_ptf$Returns, "R", pos=1)
abline(rf, max(sharpes_list))


##plotting utility function########
A <- 4
y_star <- (opt_ptf_exp_ret - rf) / (A*opt_ptf_var)
opt_c_exp_ret <- rf + y_star*(opt_ptf_exp_ret - rf)
opt_c_sqrt <- y_star*sqrt(opt_ptf_var)
u_score <-  opt_c_exp_ret - 0.5*A*opt_c_sqrt^2
sgs <- seq(0.001, 0.010, by=0.001)
rts <- sapply(sgs, function(s) {u_score + 0.5*A*s^2})
lines(sgs, rts)


##plotting optimal complete portfolio starting for the minimum variance risky portfolio

points(opt_c_sqrt, opt_c_exp_ret, col="green", cex=1, pch=16)
text(opt_c_sqrt, opt_c_exp_ret, "C", pos=1)
#######################

points(sqrt(opt_ptf_var), opt_ptf_exp_ret, col="blue", pch=16)



##### 2 SECTION TO RUN ################

#### Adjusting optimal weights for tactical asset allocation #####
margin <- 0.05

low_limits <- setNames(rep(0, n), tickers)


if (curr_rets[length(curr_rets)] > 0 && curr_rets[length(curr_rets) - 1] > 0) {
  print("RISK OFF")
  w_dict["MSFT"] <- w_dict["MSFT"] - margin
  low_limits["MSFT"] <- w_dict["MSFT"]
  w_dict["JNJ"] <- w_dict["JNJ"] + margin
  low_limits["JNJ"] <- w_dict["JNJ"]
  opt <- hydroPSO(fn = v.c, method="fips", control=list(MinMax="min", write2disk=FALSE), lower=low_limits, upper=rep(1, n))
  pos_w_star <- opt$par
  points(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, col="purple", pch=16)
  text(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, "RISK OFF")

} else if (curr_rets[length(curr_rets)] < 0 && curr_rets[length(curr_rets) - 1] < 0) {
  print("RISK ON")
  w_dict["MSFT"] <- w_dict["MSFT"] + margin
  low_limits["MSFT"] <- w_dict["MSFT"]
  w_dict["JNJ"] <- w_dict["JNJ"] - margin
  low_limits["JNJ"] <- w_dict["JNJ"]

  opt <- hydroPSO(fn = v.c, method="fips", control=list(MinMax="min", write2disk=FALSE), lower=low_limits, upper=rep(1, n))
  pos_w_star <- opt$par
  points(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, col="purple", pch=16)
  text(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, "RISK ON")
} else {print("Not enough ST information to apply a strategy") }
# ##############################################################
 




### 3 SECTION TO RUN #################
#### Adjusting optimal weights for strategical asset allocation #####

margin <- 0.02; 


scenarios <- setNames(c("disinflation-reflation", "high-inflation", "deflation"), c(1, 2, 3))
actual_scenario <- scenarios[1]


if (actual_scenario == "disinflation-reflation") {
  print("disinflation-reflation")
  w_dict["MCD"] <- w_dict["MCD"] + margin
  low_limits["MCD"] <- w_dict["MCD"]
  w_dict["JNJ"] <- w_dict["JNJ"] - margin
  low_limits["JNJ"] <- w_dict["JNJ"]
  
  opt <- hydroPSO(fn = v.c, method="fips", control=list(MinMax="min", write2disk=FALSE), lower=low_limits, upper=rep(1, n))
  pos_w_star <- opt$par
  points(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, col="purple", pch=16)
  text(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, "DR")
} else if (actual_scenario == "high-inflation")  {
  print("high-inflation")
  w_dict["MCD"] <- w_dict["MCD"] - margin
  low_limits["MCD"] <- w_dict["MCD"]
  w_dict["JNJ"] <- w_dict["JNJ"] + margin
  low_limits["JNJ"] <- w_dict["JNJ"]
  
  opt <- hydroPSO(fn = v.c, method="fips", control=list(MinMax="min", write2disk=FALSE), lower=low_limits, upper=rep(1, n))
  pos_w_star <- opt$par
  points(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, col="purple", pch=16)
  text(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, "INF")
} else if (actual_scenario == "deflation") {
  print("deflation")
  w_dict["MCD"] <- w_dict["MCD"] - margin
  low_limits["MCD"] <- w_dict["MCD"]
  w_dict["JNJ"] <- w_dict["JNJ"] - margin
  low_limits["JNJ"] <- w_dict["JNJ"]
  
  opt <- hydroPSO(fn = v.c, method="fips", control=list(MinMax="min", write2disk=FALSE), lower=low_limits, upper=rep(1, n))
  pos_w_star <- opt$par
  points(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, col="purple", pch=16)
  text(sqrt(t(opt$par)%*% S %*% opt$par), t(as.matrix(opt$par)) %*% r, "DEF")
}

###########################


##plotting rf label
text(0.0002, rf, "rf")





