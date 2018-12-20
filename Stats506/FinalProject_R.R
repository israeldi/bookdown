# Final Project STATS 506
# Monte Carlo Simulation of Portfolio Stock Returns
#
# This is a tutorial code for running a simple monte carlo simulation to
# estimate the value of our Portfolio one period ahead.
#
# Apple, Google, and Facebook adjusted close Prices from November 14th, 2017 - 
# November 14th, 2018. Data comes from Yahoo Finance.
#
# Author: Shuoran Li (shuoranl@umich.edu)
# Due Date: December 7, 2018

# Load library
library(data.table)
library(ggplot2)
# Load data 
stock_Data = fread('./Group21_ProjectData.csv')
stock_Price = as.matrix( stock_Data[ , 2:4] )

mc_rep = 1000 # Number of Monte Carlo Simulations
training_days = 30 

# This function returns the first differences of a t x q matrix of data
returns = function(Y){
  len = nrow(Y)
  yDif = Y[2:len, ] / Y[1:len-1, ] - 1
}

# Get the Stock Returns
stock_Returns = returns(stock_Price)

# Suppose we invest our money evenly among all three assets 
# We use today's Price 11/14/2018 to find the number of shares each stock 
# that we buy
portfolio_Weights = t(as.matrix(rep(1/ncol(stock_Returns), 
                                    ncol(stock_Returns))))


# Get the Variance Covariance Matrix of Stock Returns
coVarMat = cov(stock_Returns)
miu = colMeans(stock_Returns)
# Extend the vector to a matrix
Miu = matrix(rep(miu, training_days), nrow = 3)

# Initializing simulated 30 day portfolio returns
portfolio_Returns_30_m = matrix(0, training_days, mc_rep)

set.seed(200)
for (i in 1:mc_rep) {
  Z = matrix ( rnorm( dim(stock_Returns)[2] * training_days ), 
               ncol = training_days )
  # Lower Triangular Matrix from our Choleski Factorization
  L = t( chol(coVarMat) )
  # Calculate stock returns for each day
  daily_Returns = Miu + L %*% Z  
  # Calculate portfolio returns for 30 days
  portfolio_Returns_30 = cumprod( portfolio_Weights %*% daily_Returns + 1 )
  # Add it to the monte-carlo matrix
  portfolio_Returns_30_m[,i] = portfolio_Returns_30;
}

# Visualising result
x_axis = rep(1:training_days, mc_rep)
y_axis = as.vector(portfolio_Returns_30_m - 1)
plot_data = data.frame(x_axis, y_axis)
ggplot(data = plot_data, aes(x = x_axis, y = y_axis)) + 
  geom_path(aes(col = 'red'), size = 0.1) +
  xlab('Days') + ylab('Portfolio Returns') + 
  ggtitle('Simulated Portfolio Returns in 30 days')


# Porfolio Returns statistics at the 30th day.

Avg_Portfolio_Returns = mean(portfolio_Returns_30_m[training_days,] - 1)
SD_Portfolio_Returns = sd(portfolio_Returns_30_m[training_days,] - 1)
Median_Portfolio_Returns = median(portfolio_Returns_30_m[training_days,] - 1)

# Construct a 95% Confidential Interval for average returns

Avg_CI = quantile(portfolio_Returns_30_m[training_days,] - 1, c(0.025, 0.975))





