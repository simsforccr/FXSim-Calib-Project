
# coding: utf-8

# In[1]:


import os
os.chdir('/home/malek/Code/FXSim-Calib-Project/exercises/')
os.listdir()


# In[43]:


import pandas as pd
rates = pd.read_csv('eurofxref-hist.csv',usecols=['Date','USD','GBP','JPY','CHF','AUD','CAD'])
print(rates.head())
rates.Date = pd.to_datetime(rates.Date)
rates[rates.Date > '2018-01-01'].head()
rates = rates.set_index('Date')


# In[49]:


rates['EUR'] = 1
rates = rates[rates.index < '2016-01-01']
rates = rates[rates.index > '2012-12-31']
rates = rates.div(rates.GBP,axis=0)
rates = rates.drop('GBP',1)


# In[74]:


#we calculate log returns from the past 3 years
logreturns = np.log(rates.USD) - np.log(rates.shift(-1).USD)
#the standard deviation of our simulated log returns will be the standard deviations of the historical data
std = logreturns.std()
nsims = 100
sim_dates = pd.bdate_range('2015-12-31','2016-12-31')
import numpy as np

#we generate normal random numbers, which will be the log returns of our simulated prices.
rnds = np.random.normal(0,std,[len(sim_dates),nsims])
rnorm = pd.DataFrame(rnds,columns=range(nsims),index=sim_dates)


# In[87]:


sims = pd.DataFrame(columns = range(nsims),index = sim_dates)

#our simulated data will start with the last data point from our historical data
start_date = sims.index[0]
sims.loc[start_date] = rates.USD[start_date]


# In[88]:


for i in range(1, len(sim_dates)):
    sims.loc[sim_dates[i]] = np.exp(rnorm.loc[sim_dates[i]])*sims.loc[sim_dates[i-1]]

print(sims.head())
get_ipython().run_line_magic('matplotlib', 'inline')
sims.plot(legend = False)

