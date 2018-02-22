Please set you directory in the main program eurusd.R

You can either run it on one of the dataset by setting the global NEW to a logical variable or run it from scratch by setting NEW to a string. Be warned that you need a connection to Interactive brokers should you want to download new data.

The program can have difficulties in running on Windows. My current setup was on bash on Windows. i7 core with 8GB of RAM.

Although it uses parallel programming, the program can be very slow due to optimization of algorithms.

The program consists of:
- Downloading all G7 currencies from Interactive brokers
- Fitting arima garch on each of the currencies
- Fitting arima monte carlo on currencies
- Include different variables such as lags and past returns
- Train extraTrees and Random Forest to the currencies selected in GLOBALS (Run on all by default - 21 currencies pair).

Compute the total performance of the strategy.
