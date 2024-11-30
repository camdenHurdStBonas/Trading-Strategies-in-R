# Clear the workspace
rm(list = ls())

# Set working directory
path <- "/Users/camdenhurd/desktop/Desktop/Projects/Robinhood Algo Trading/robinhood-api-trading"
setwd(path)

# Define the risk-free rate (annualized)
risk_free_rate <- 0.03
# Define transaction cost (e.g., 0.005 or 0.01 as per your suggestion)
transaction_cost <- 0.005
# Define reduced ranges for the parameters to speed up computations
nFast_values <- seq(20, 40, 1)    
nSlow_values <- seq(40, 80, 1)     
nSig_values <- seq(7, 10, 1)
vwap_window_values <- seq(10, 80, 1)
tema_window_values <- seq(10, 80, 1)

# Load historical data
asset_name <- "BTC-USD"
start_date <- "2017-01-01"

# Run the strategies and get their signals
source("MACD.R")
macd <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost,  nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values)
source("VWAP.R")
vwap <- vwap_strategy(asset_name = asset_name, start_date = start_date, vwap_window_values = vwap_window_values, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost)
source("MVCD.R")
mvcd <- mvcd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost,  nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values)
source("TEMA.R")
tema <- tema_strategy(asset_name = asset_name, start_date = start_date, tema_window_values = tema_window_values, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost)

# Combine signals into a matrix
strategy_returns <- cbind(macd$best_strategy$returns,vwap$returns,tema$returns,mvcd$best_strategy$returns)
strategy_names <- c("MACD","VWAP","TEMA","MVCD")
colnames(strategy_returns) <- strategy_names
# Combine signals into a matrix
strategy_signals <- cbind(macd$best_strategy$signal, vwap$signal, tema$signal,mvcd$best_strategy$signal)
strategy_names <- c("MACD","VWAP","TEMA","MVCD")
colnames(strategy_signals) <- strategy_names

# Define multiple thresholds for buying and selling
buy_thresholds <- seq(0.05, 0.5, by = 0.05)
sell_thresholds <- seq(-.5, -0.05, by = 0.05)
stop_loss_percentages <- seq(0.15, 0.25,by = 0.05)
take_profit_percentages <-  seq(8.0, 10.0,by = 0.5)
risk_percentages <- seq(0.05, 0.10, by = 0.01)

source("optimize_trading_indicator_weights.R")
optimal_weights <- optimize_weights(asset_name,start_date,strategy_returns,strategy_signals, stock_returns, risk_free_rate,transaction_cost, strategy_names, buy_thresholds, sell_thresholds, stop_loss_percentages,take_profit_percentages,risk_percentages)

cryptos <- c(
  "BTC-USD",  # Bitcoin
  "ETH-USD",  # Ethereum
  "LTC-USD",  # Litecoin
  "BCH-USD",  # Bitcoin Cash
  "DOGE-USD", # Dogecoin
  "ETC-USD",  # Ethereum Classic
  "BSV-USD",  # Bitcoin SV
  "COMP5692-USD", # Compound
  "MATIC-USD",# Polygon
  "SOL-USD",  # Solana
  "SHIB-USD", # Shiba Inu
  "LINK-USD", # Chainlink
  "ADA-USD",  # Cardano
  "XRP-USD",  # XRP
  "PEPE24478-USD", # Pepe
  "XLM-USD",  # Stellar
  "AVAX-USD", # Avalanche
  "UNI7083-USD",  # Uniswap
  "XTZ-USD",   # Tezos
  "WIF-USD"    # Dogwifhat
)
start_date <- "2024-06-01"
source("market_performance.R")
market_performance <- market_performance(cryptos, start_date, rolling_window = 30) 





  










optimized_params <- list(
  optimal_weights = optimal_weights$OptimalWeights,
  nFast = macd$nFast,                # MACD fast period
  nSlow = macd$nSlow,                # MACD slow period
  nSig = macd$nSig,                  # MACD signal period
  vwap_window = vwap$BestVWAPWindow, # VWAP best window
  tema_window = tema$BestTEMAWindow, # TEMA best window
  mvcd_nFast = mvcd$nFast,           # MVCD fast period
  mvcd_nSlow = mvcd$nSlow,           # MVCD slow period
  mvcd_nSig = mvcd$nSig,             # MVCD signal period
  buy_threshold = optimal_weights$best_buy_threshold,
  sell_threshold = optimal_weights$best_sell_threshold,
  stop_loss = optimal_weights$BestTakeProfit,
  take_profit = optimal_weights$BestStopLoss 
)


source("validate_parameters.R")
validation_results <- validate_parameters(
  asset_name = asset_name,
  start_date = start_date,
  optimized_params = optimized_params,
  risk_free_rate = risk_free_rate,
  transaction_cost = transaction_cost
)




source("MACD_Advanced.R")
macd_advanced_strategy(asset_name, start_date, risk_free_rate, 0.05, transaction_cost, nFast_values, nSlow_values, nSig_values, short = 0)











asset_name <- "ETH-USD"
start_date <- "2019-10-01"
strategy_2 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)
macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)
rsi_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, rsi_periods = rsi_periods, upper_thresholds = upper_thresholds, lower_thresholds = lower_thresholds, short = short)

asset_name <- "SOL-USD"
start_date <- "2020-09-01"
strategy_3 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "DOGE-USD"
start_date <- "2019-10-01"
strategy_4 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "XRP-USD"
start_date <- "2019-10-01"
strategy_5 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "ADA-USD"
start_date <- "2019-10-01"
strategy_6 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "SHIB-USD"
start_date <- "2021-07-01"
strategy_7 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "AVAX-USD"
start_date <- "2020-10-01"
strategy_8 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "LINK-USD"
start_date <- "2019-10-01"
strategy_9 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "BCH-USD"
start_date <- "2019-10-01"
strategy_10 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "UNI7083-USD"
start_date <- "2020-10-01"
strategy_11 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "PEPE24478-USD"
start_date <- "2023-06-01"
strategy_12 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "LTC-USD"
start_date <- "2019-10-01"
strategy_13 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "XLM-USD"
start_date <- "2019-10-01"
strategy_14 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "ETC-USD"
start_date <- "2019-10-01"
strategy_15 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "AAVE-USD"
start_date <- "2020-11-01"
strategy_16 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "XTZ-USD"
start_date <- "2019-10-01"
strategy_17 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)

asset_name <- "COMP5692-USD"
start_date <- "2020-07-01"
strategy_18 <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)



result <- rsi_strategy(asset_name = asset_name, start_date = "2022-01-01", risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowing_cost = borrowing_cost, rsi_periods = rsi_periods, upper_thresholds = upper_thresholds, lower_thresholds = lower_thresholds, short = short)
















# Define the parameters
asset_name <- "AAPL"
start_date <- "2018-01-01"
risk_free_rate <- 0.03   # Annualized risk-free rate
transaction_cost <- 0.005 # Transaction cost per trade
borrowing_cost <- .05 # Annualized cost of borrowing
nFast_values <- seq(8, 12, 1)
nSlow_values <- seq(20, 30, 1)
nSig_values <- seq(7, 9, 1)
short <- -1  # 0 = only buy, -1 = buy and short, 1 = only short

# Run the MACD strategy
strategy <- macd_strategy(asset_name = asset_name, start_date = start_date, risk_free_rate = risk_free_rate, transaction_cost = transaction_cost, borrowinf_cost <- borrowing_cost, nFast_values = nFast_values, nSlow_values = nSlow_values, nSig_values = nSig_values, short = short)















# Load necessary libraries
library(quantmod)
library(PerformanceAnalytics)
library(TTR)
library(ggplot2)
library(parallel)
library(plotly)

asset_name <- "BTC-USD"

start_date <- "2019-01-01"

# Define the risk-free rate (annualized)
risk_free_rate <- 0.03

# Define transaction cost (e.g., 0.005 or 0.01 as per your suggestion)
transaction_cost <- 0.005

# Define reduced ranges for the parameters to speed up computations
nFast_values <- seq(1, 30, 1)    
nSlow_values <- seq(1, 60, 1)     
nSig_values <- seq(3, 5, 1)




# Download stock data
getSymbols(asset_name, from = start_date, to = Sys.Date(), src = "yahoo")

# Get closing prices using the asset name
stock_data <- na.omit(Cl(get(asset_name)))  # Closing price of stock
stock_returns <- dailyReturn(stock_data)

# Function to evaluate MACD with given parameters
evaluate_macd <- function(nFast, nSlow, nSig, rf_rate, transaction_cost) {
  # Calculate MACD
  macd_values <- MACD(stock_data, nFast = nFast, nSlow = nSlow, nSig = nSig, maType = "EMA")
  
  # Generate trading signals
  signal <- Lag(ifelse(macd_values$macd > macd_values$signal, 1, -1))
  signal[is.na(signal)] <- 0
  signal <- na.locf(signal)
  
  # Add transaction cost whenever there's a change in position
  signal_shift <- Lag(signal)  # Previous day's signal to detect changes
  signal_shift[is.na(signal_shift)] <- 0
  trade_occurred <- signal != signal_shift  # Identify trade days
  
  # Calculate returns from the trading strategy
  strategy_returns <- stock_returns * Lag(signal)
  strategy_returns <- strategy_returns - transaction_cost * trade_occurred  # Subtract transaction cost on trade days
  
  # Calculate Sharpe Ratio for the strategy (subtract risk-free rate)
  excess_returns <- log(1 + strategy_returns) - (rf_rate / 252)  # Adjusting for daily returns
  sharpe_ratio <- SharpeRatio.annualized(excess_returns, scale = 252)
  
  return(sharpe_ratio)
}

# Data frame to store results
results <- expand.grid(nFast = nFast_values, nSlow = nSlow_values, nSig = nSig_values)
results$SharpeRatio <- NA

# Remove rows where nFast is greater than or equal to nSlow and nSig is greater than nFast
results <- subset(results, nFast < nSlow & nSig <= nFast)

# Use mclapply for parallel processing
results$SharpeRatio <- mclapply(1:nrow(results), function(i) {
  nFast <- results$nFast[i]
  nSlow <- results$nSlow[i]
  nSig <- results$nSig[i]
  evaluate_macd(nFast, nSlow, nSig, risk_free_rate, transaction_cost)
}, mc.cores = detectCores() - 1)  # Use available cores minus one

# Ensure SharpeRatio is numeric
results$SharpeRatio <- as.numeric(results$SharpeRatio)

# Step 5: 3D Visualization of MACD Parameter Combinations

# Create a 3D plot with Sharpe Ratio on the Z-axis
plot_ly(data = results, 
        x = ~nFast, 
        y = ~nSlow, 
        z = ~SharpeRatio, 
        type = 'scatter3d', 
        mode = 'markers', 
        marker = list(size = 3), 
        color = ~SharpeRatio,  # Color based on Sharpe Ratio
        colors = colorRamp(c("red", "green"))) %>%
  layout(title = "3D Plot of Sharpe Ratio by nFast and nSlow",
         scene = list(xaxis = list(title = "nFast"),
                      yaxis = list(title = "nSlow"),
                      zaxis = list(title = "Sharpe Ratio")),
         legend = list(title = list(text = "Sharpe Ratio")))

# Step 6: Recalculate Best MACD Strategy and Comparison

# A. Identify the best parameters
best_result <- results[which.max(results$SharpeRatio), ]
cat("Best Parameters:\n")
print(best_result)

# Best parameters
best_nFast <- best_result$nFast
best_nSlow <- best_result$nSlow
best_nSig <- best_result$nSig

# B. Calculate MACD
macd_values <- MACD(stock_data, nFast = best_nFast, nSlow = best_nSlow, nSig = best_nSig, maType = "EMA")

# C. Create a trading strategy based on MACD crossovers
# Buy when MACD crosses above the Signal line, sell when it crosses below
signal <- Lag(ifelse(macd_values$macd > macd_values$signal, 1, -1))  # Generate trading signals
signal[is.na(signal)] <- 0  # Replace NA values
signal <- na.locf(signal)  # Carry last signal forward

# D. Add transaction cost for signal changes
signal_shift <- Lag(signal)  # Previous signal
signal_shift[is.na(signal_shift)] <- 0
trade_occurred <- signal != signal_shift

# D. Calculate returns from the trading strategy
strategy_returns <- stock_returns * Lag(signal)  # Strategy returns based on MACD signals
strategy_returns <- strategy_returns - transaction_cost * trade_occurred

# E. Buy-and-hold strategy returns
buy_hold_returns <- stock_returns  # Buy and hold returns are just holding stock

# F. Plot MACD and Signal with the Price
chartSeries(stock_data, TA=paste("addMACD(", best_nFast, ", ", best_nSlow, ", ", best_nSig, ")"))

# G. Compare MACD Strategy vs Buy-and-Hold Strategy

# Combine the returns from both strategies into a single object
comparison_returns <- cbind(buy_hold_returns, strategy_returns)
colnames(comparison_returns) <- c("Buy and Hold", "MACD Strategy")

# Plot performance of both strategies
charts.PerformanceSummary(comparison_returns, 
                          main = "MACD Strategy vs. Buy and Hold", 
                          legend.loc = "topleft", 
                          colorset = c("black", "blue"))

# H. Output basic performance metrics for both strategies
annualized_returns <- table.AnnualizedReturns(comparison_returns, scale = 252)  # 252 is used for daily returns

# Extract the metrics for Buy and Hold and MACD Strategy
annualized_return_bh <- annualized_returns["Annualized Return", "Buy and Hold"]
annualized_return_macd <- annualized_returns["Annualized Return", "MACD Strategy"]
annualized_std_dev_bh <- annualized_returns["Annualized Std Dev", "Buy and Hold"]
annualized_std_dev_macd <- annualized_returns["Annualized Std Dev", "MACD Strategy"]

# Subtract risk-free rate for Sharpe Ratio calculation
excess_buy_hold <- log(1+comparison_returns[, "Buy and Hold"]) - (risk_free_rate / 252)
excess_macd <- log(1+comparison_returns[, "MACD Strategy"]) - (risk_free_rate / 252)

# Calculate Sharpe Ratios
sharpe_bh <- SharpeRatio.annualized(excess_buy_hold, scale = 252)
sharpe_macd <- SharpeRatio.annualized(excess_macd, scale = 252)

# Calculate Max Drawdown for each strategy
max_drawdown_bh <- maxDrawdown(comparison_returns[,"Buy and Hold"])
max_drawdown_macd <- maxDrawdown(comparison_returns[,"MACD Strategy"])

# Combine all metrics into a data frame
performance_metrics <- data.frame(
  Metric = c("Annualized Return", "Annualized Std Dev", "Sharpe Ratio", "Max Drawdown"),
  Buy_and_Hold = c(annualized_return_bh, annualized_std_dev_bh, sharpe_bh, max_drawdown_bh),
  MACD_Strategy = c(annualized_return_macd, annualized_std_dev_macd, sharpe_macd, max_drawdown_macd)
)

# I. Print the combined performance metrics table
print(performance_metrics)
