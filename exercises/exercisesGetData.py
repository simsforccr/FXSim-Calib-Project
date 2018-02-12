import pandas as pd
from scipy.linalg import cholesky
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

def ex1_results():
    print rates.head()

#ex1_results()


#we will work mostly with log returns of the FX pairs. 
logreturns = np.log(rates)-np.log(rates.shift(-1))

def ex2_results():
    print logreturns.corr()
    
#ex2_results()


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

def ex3_results():
    print sims.head()
    sims.plot(legend= False)


#ex3_results()

class FXForward:
    def __init__(self,pay_ccy='USD',rec_ccy='GBP',pay_amt=1500000,rec_amt=1000000,pay_date=pd.to_datetime('2016-10-10'),rec_date=pd.to_datetime('2016-10-10')):
        self.pay_ccy = pay_ccy
        self.rec_ccy = rec_ccy
        self.rec_amt=rec_amt
        self.pay_amt = pay_amt
        self.pay_date = pd.to_datetime(pay_date)
        self.rec_date = pd.to_datetime(rec_date)
        
    def __repr__(self):
        result = ''
        result += 'pay: '+ self.pay_ccy + ' ' + str(self.pay_amt) + '\n'
        result += 'rec: '+ self.rec_ccy + ' ' + str(self.rec_amt) + '\n'
        result += str(self.pay_date) + '\n'
        
        


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
        #we'll have to change this when we want to be able to h trades in more than one currency
       
    
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

def ex4_results():
    fxforward = FXForward()
    print fxforward

def ex5_results():
    prices = fxforward.price_from_sims(sims)
    print prices.head()
    prices.plot()

def exposure_profile(trades,sims):
    mtm = pd.DataFrame(0,index = sims.index,columns = sims.columns)
    for trade in trades:
        mtm = mtm + trade.price_from_sims(sims)
    
    return mtm

def epe(sims):
    #takes as input a date-indexed DataFrame, and returns a date-indexed Series with average positive values. 
    return sims[sims>0].fillna(0).mean(axis=1)
    
def pfe(sims,percentile):
    return sims[sims>0].fillna(0).quantile(percentile/100.0, axis=1)
    
def make_random_usdgbp_trades(num_trades):

    trades = []
    for i in range(num_trades):
        date = np.random.choice(sim_dates)
        amount = np.random.randint(1,1000)*1000
        pay = np.random.choice(['USD','GBP'])
        rec = 'USD' if pay == 'GBP' else 'GBP'
        rate = np.random.uniform(1.45,1.5)
        usd_amt = rate*amount
        pay_amt = usd_amt if pay == 'USD' else amount
        rec_amt = usd_amt if rec == 'USD' else amount
        trade = FXForward(pay_ccy =pay,rec_ccy=rec,rec_amt = rec_amt,pay_amt=pay_amt,pay_date = date,rec_date = date)
        trades.append(trade)
    return trades

def ex9_results():
    trades = make_random_usdgbp_trades(100)
    sim_values = exposure_profile(trades,sims)

    sim_values.plot(legend=False)
    epe(sim_values).plot()
    pfe(sim_values,98).plot()


#We now return to the case of mutiple currencies. 
corr = logreturns.corr()
std = logreturns.std()
chol = cholesky(corr)
currencies = corr.columns
multi = pd.MultiIndex.from_product([currencies,sim_dates],names=['currency','date'])
sims = pd.DataFrame(index = multi, columns = range(nsims))



