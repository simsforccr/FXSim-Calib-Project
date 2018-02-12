import pandas as pd
import numpy as np
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

nsims = 100
sim_dates = pd.bdate_range('2015-12-31','2016-12-31')

#we generate normal random numbers for each of the simulated values
rnds = np.random.normal(0,std,[len(sim_dates),nsims])
df_rnorm = pd.DataFrame(rnds, columns = range(nsims),index=sim_dates)

sims = pd.DataFrame(columns= range(nsims),index=sim_dates)

last_date = rates.USD.index[0]
sims.loc[last_date] = rates.USD[last_date]
for i in range(1,len(sim_dates)):
    sims.loc[sim_dates[i]] = np.exp(df_rnorm.loc[sim_dates[i]])*sims.loc[sim_dates[i-1]]

    
sims.plot()


class FXForward:
    def __init__(self,pay_ccy='USD',rec_ccy='GBP',pay_amt=1500000,rec_amt=1000000,pay_date=pd.to_datetime('2016-10-10'),rec_date=pd.to_datetime('2016-10-10')):
        self.pay_ccy = pay_ccy
        self.rec_ccy = rec_ccy
        self.rec_amt=rec_amt
        self.pay_amt = pay_amt
        self.pay_date = pd.to_datetime(pay_date)
        self.rec_date = pd.to_datetime(rec_date)
        


    def price(self,pay_rate,rec_rate):
        if self.pay_ccy == 'GBP':
            self.pay_value = self.pay_amt
        
        else:
            self.pay_value = self.pay_amt/pay_rate
            
        if self.rec_ccy == 'GBP':
            self.rec_value = self.rec_amt
            
        else:
            self.rec_value = self.rec_amt/rec_rate
            
        self.value = self.rec_value - self.pay_value
        
        return self.value
        
    def price_from_sims(self,sims):
        #for now, sims is a nsims x sim_dates DataFrame, which will give exchange rates for the currency pair we're interested in
        #we'll have to change this when we want to be able to price trades in more than one currency
       
    
        self.pay_values = pd.DataFrame(self.pay_amt,index = sims.index,columns = sims.columns)
        self.rec_values = pd.DataFrame(self.rec_amt,index = sims.index,columns = sims.columns)
        
        #we will use the simulated FX rate at the tenor date to represent our simulation of the forward rate at that tenor. 
        
        if self.pay_ccy == 'GBP':
            pass 
        else:
            self.pay_values = self.pay_values.div(sims)
           
        self.pay_values[self.pay_values.index > self.pay_date] = 0
        
        if self.rec_ccy == 'GBP':
            pass

        else:
            self.rec_values = self.rec_values.div(sims)
        
        self.rec_values[self.rec_values.index > self.rec_date] = 0
        
        self.values = self.rec_values-self.pay_values
        return self.values

fxforward = FXForward('abc')
print fxforward.price_from_sims(sims).head()