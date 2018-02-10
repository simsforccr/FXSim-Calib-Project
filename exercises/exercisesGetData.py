mport pandas as pd
import numpy as np
from ggplot import *
#get data and load relevant columns into a pandas DataFrame
filename = 'eurofxref-hist.csv'
allrates = pd.read_csv(filename,usecols=['Date','USD','GBP','JPY','CAD','CHF','AUD'])

#we could do this in the call to read_csv, but we want to the DataFrame to know the dates are dates.
allrates['Date'] = pd.to_datetime(allrates['Date'])

#we will use the 'Date' column as an index.
allrates = allrates.set_index('Date')

#we want to use GBP as our base currency, so we add a 'EUR' column, and divide everthing by the GBP rate.
allrates['EUR'] = 1
allrates = allrates.div(allrates.GBP,axis=0)

#this will throw off some later calculations if we keep it (as it never changes), so we drop it from the DataFrame.
allrates = allrates.drop('GBP',1)


#We only want three years of history.
rates = allrates[allrates.index <= '2015-12-31']
rates = rates[rates.index > '2012-12-31']

#we will work mostly with log returns of the FX pairs. 
logreturns = np.log(rates)-np.log(rates.shift(-1))


corr = logreturns.corr()

#for the first day, we will focus on the USD/GBP rate.
usd = logreturns.USD
std = usd.std()

nsims = 1000
sim_dates = pd.bdate_range('2015-12-31','2016-12-31')

#we generate normal random numbers for each of the simulated values
rnds = np.random.normal(0,std,[nsims,len(sim_dates)])
df_rnorm = pd.DataFrame(rnds, index = range(nsims),columns=sim_dates)

sims = pd.DataFrame(index = range(nsims),columns=sim_dates)

last_date = rates.USD.index[0]
sims[last_date] = rates.USD[last_date]
for i in range(1,len(sim_dates)):
    sims[sim_dates[i]] = np.exp(df_rnorm[sim_dates[i]])*sims[sim_dates[i-1]]

    
ggplot(sims)


