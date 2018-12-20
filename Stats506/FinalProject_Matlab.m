% Final Project STATS 506
% Monte Carlo Simulation of Portfolio Stock Returns
%
% This is a tutorial code for running a simple monte carlo simulation to
% estimate the value of our Portfolio one period ahead.
%
% Apple, Google, and Facebook adjusted close Prices from November 14th, 2017 - 
% November 14th, 2018. Data comes from Yahoo Finance.
%
% Author: Israel Diego (israeldi@umich.edu)
% Due Date: December 7, 2018

clear
clc

% addpath '/Users/israeldiego/Documents/GitHub/506_Project'

% Set Seed of Random Number Generator 
rng('default');
rng(1);

mc_rep = 1000; % Number of Monte Carlo Simulations
initInvestment = 100000; % $100,000 Initial Investment Portfolio at t = 0
numTradingDays = 30;

% Load Stock Price Data
stockData = readtable('Group21_ProjectData.csv');
stockPrices = table2array(stockData(:, 2:end));

% Calculate our Daily Returns
stockReturns = returns(stockPrices);

% Suppose we invest our money evenly among all three assets 
% We use today's Price 11/14/2018 to find the number of shares each stock
% that we buy
portfolioWeights = (1/3) * ones(1, size(stockPrices,2));

% Get the Variance Covariance Matrix of our Stock Returns
coVarMat = cov(stockReturns);

% Average returns of each asset 
mu = transpose(mean(stockReturns));
mu = repmat(mu, 1, numTradingDays);

% Initializing our simulated 30 day returns
portfolio30DayReturn_m = zeros(numTradingDays, mc_rep);
for i = 1:mc_rep
    % Randomly generated numbers from N(0,1) distribution 
    Z = randn(size(stockReturns,2), numTradingDays);

    % Lower Triangular Matrix from Choleski Factorization
    L = chol(coVarMat, 'lower');

    % Calculate daily returns for 30 days
    dailyReturns = mu + (L * Z);

    % Portfolio Returns
    thirtyDayReturn = transpose(cumprod(portfolioWeights * dailyReturns + 1));
    
    % Add return to the set of all 30-day portfolio returns
    portfolio30DayReturn_m(:,i) = thirtyDayReturn;
end

plot(portfolio30DayReturn_m - 1, 'LineWidth',0.5, 'Color',[0,0.7,0.9, 0.2])
title('Simulated Portfolio Returns in 30 days', 'fontsize', 16)
xlabel('Days','fontsize',16)
ylabel('Portfolio Returns','fontsize',16)

% Save Plot
saveas(gcf,[pwd '/matlab_pics/ThirtyDay_Returns_Matlab_MC.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate some statistics for our simulated portfolio returns
averagePortfolioReturns = mean(portfolio30DayReturn_m(end,:) - 1);
stdDevPortfolioReturns = std(portfolio30DayReturn_m(end,:) - 1);
medianPortfolioReturns = median(portfolio30DayReturn_m(end,:) - 1);

% Construct a 95% Confidential Interval for average returns

average_CI = quantile(portfolio30DayReturn_m(end,:) - 1, [0.025, 0.975]);

% Display Results
display(['Statistics of Simulated 30-day Portfolio Returns']);
fprintf('\n')
display(['------------------------------------------------']);
fprintf('\n')
fprintf('\n')
display(['Mean 30-day Returns: ' num2str(averagePortfolioReturns)]);
fprintf('\n')
display(['Standard Deviation 30-day Returns: ' ...
    num2str(stdDevPortfolioReturns)]);
fprintf('\n')
display(['Median 30-day Returns: ' num2str(medianPortfolioReturns)]);
fprintf('\n')
display(['95% CI 30-day Returns: ' num2str(average_CI)]);
fprintf('\n')

% This function returns the first differences of a t x q matrix of data
function [yDif] = returns(y)
    yDif = y(2:end, :) ./ y(1:end-1, :) - 1;
end































