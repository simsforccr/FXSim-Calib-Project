# -*- coding: utf-8 -*-
"""
Soluce pour FX calib, sim and pricing

M.Jawad - 04/01/2018

"""
import datetime
import pandas as pd
import numpy as np
import scipy as sp
from scipy import linalg

from matplotlib import pyplot as plt

Path = "C:\\Users\\Malek\\Documents\\Python Projects\\FXSim-Calib-Project\\"
startDate = datetime.date(2015,1,2)
endDate = datetime.date(2015,12,31)

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + datetime.timedelta(n)

NbSims = 1000
SimLength = 365
Percentile = 99

def FXSim(fxspot, vol, rdm, time):
    """
    FX simulation for n days with x vol using the random numbers provided over the "time" horizon
    """
    # create an array of lenght NbSims x "time" horizon
    sim = np.zeros((NbSims,len(time)))
    
    # simulate a random walk FX with no IR discounting
    for i in range(0,NbSims):
        for j in range(0,len(time)):
            sim[i,j] = fxspot+time[j]*vol*rdm[i,j]
    
    # return simulation
    return sim

def SimulateFXRates(FilePath, SD, ED, NbSimulations, SimHorizon):
    """
    main function to simulate FX rates for each day
    input is csv file path, startdate, enddate, Number of simulations and simHorizon (nb of days)
    
    output 1 is a 4 dimension numpy array shaped by n days, currencyIndex, Nb simulations,  simHorizon (nb of days)
    output 2 is the currency list as reference
    output 3 is a dataframe of dates for reference for valuation
    """
    
    # load historical timeseries
    dateparse = lambda x: pd.datetime.strptime(x, '%d/%m/%Y')
    df = pd.read_csv(FilePath, parse_dates=['DATE'], date_parser=dateparse)

    # setup sim variables
    SimCalibration = (ED - SD).days
    time_array = np.linspace(0,1,SimCalibration)
    CcyList = list(df)[1:]
    
    # generate log returns dateframe and rolling volatility over the length of the calibration period
    df_LogR = np.log(df.loc[:,CcyList]) - np.log(df.loc[:,CcyList].shift(1))
    df_Vol = df_LogR.rolling(SimCalibration, SimCalibration).std()*sp.sqrt(SimCalibration)
    df_Vol = pd.concat([df.loc[:,['DATE']],df_Vol], axis=1)

    dateRange = [df_Vol.index[df_Vol['DATE'] == SD].tolist()[0], df_Vol.index[df_Vol['DATE'] == ED].tolist()[0]]
    SimArray = np.zeros((dateRange[1]-dateRange[0]+1,len(CcyList),NbSimulations,len(time_array)))
    # generate daily FX simulations over 1 year
    for n in range(dateRange[0],dateRange[1]+1):
        
        # generate the correlation matrix from the log return dataframe for each day
        df_Corr = df_LogR.loc[(n-dateRange[0]):n,CcyList].corr(method='pearson')
        
        # generate for every n day correlated random numbers til the SimHorizon
        CorrRdm = np.zeros((len(CcyList),NbSimulations,SimHorizon))
        for j in range(0, SimHorizon):
            CorrRdm[:,:,j] = np.dot(linalg.cholesky(df_Corr) , np.random.normal(0, 1, size=(len(CcyList),NbSimulations)))
            
        # generate FX simulations for every n day using the Correlated Random numbers generate before
        for ccy in CcyList:
            ccyI = CcyList.index(ccy)
            SimArray[n-dateRange[0],ccyI,:,:] = FXSim(df.loc[n, ccy], df_Vol.loc[n,ccy], CorrRdm[ccyI,:,:], time_array)           

    return [SimArray, CcyList, df.loc[:,['DATE']]]


# define a trade class
class FXfwdTrade:
    def __init__(self,TradeStart,MatDate,RecNot,RecCcy,PayNot,PayCcy):
        self.TradeStartDate = TradeStart
        self.maturityDate = MatDate   
        self.RecLegNotional = RecNot
        self.RecLegCcy = RecCcy
        self.RecLegCFDate = self.maturityDate
        self.PayLegNotional = PayNot
        self.PayLegCcy = PayCcy
        self.PayLegCFDate = self.TradeStartDate  

    def GenerateMTF(self, BatchDate, Dates, CcyList, Sims):
        if np.datetime64(startDate) in Dates.values and np.datetime64(BatchDate) in Dates.values:
            if BatchDate > self.maturityDate:
                self.MTF = np.zeroes((NbSims,len(CcyList)))
                            
            BatchStartIndex = Dates.index[Dates['DATE'] == startDate].tolist()[0]
            BatchIndex = Dates.index[Dates['DATE'] == BatchDate].tolist()[0] - BatchStartIndex
            MaturityIndex = (self.maturityDate - BatchDate).days
            
            # get the receive leg
            if self.RecLegCcy == 'GBP':
                RecGBPNot = self.RecLegNotional * np.ones(np.shape(Sims[BatchIndex,0,:,:MaturityIndex]))
                #* np.exp(-MaturityIndex/365*DF[BatchIndex, RecCcyIndex,:,MaturityIndex])
            else:
                RecCcyIndex = CcyList.index(self.RecLegCcy)
                RecGBPNot = self.RecLegNotional/Sims[BatchIndex,RecCcyIndex,:,:MaturityIndex] #* np.exp(-MaturityIndex/365*DF[BatchIndex, RecCcyIndex,:,:MaturityIndex])
            
            # get the pay leg
            if self.PayLegCcy == 'GBP':
                PayGBPNot = self.PayLegNotional * np.ones(np.shape(Sims[BatchIndex,0,:,:MaturityIndex]))
                #* np.exp(-MaturityIndex/365*DF[BatchIndex, RecPayIndex,:,MaturityIndex])
            else:
                PayCcyIndex = CcyList.index(self.PayLegCcy)
                PayGBPNot = self.PayLegNotional/Sims[BatchIndex,PayCcyIndex,:,:MaturityIndex] #* np.exp(-MaturityIndex/365*DF[BatchIndex, PayCcyIndex,:,:MaturityIndex]
            # price the forward or spot from RecGBPNot and PayGBPNot (both are numpy array 1000 x days to maturity)
            self.MTF = RecGBPNot - PayGBPNot
        else:
            if np.datetime64(startDate) not in Dates.values:
                print("startDate not a trading day")
            if np.datetime64(BatchDate) not in Dates.values:
                print("BatchDate not a trading day")

    def MTM(self):
        if self.MTF is not None:
            return np.average(self.MTF[:,-1])

    def EE(self):
        if self.MTF is not None:
            E = self.MTF[:,:]
            E[E < 0] = 0
            return np.mean(E, axis=0)
        
    def PFE(self,Percent):
        if self.MTF is not None:
            return np.percentile(self.MTF[:,:],Percent,axis=0,interpolation='nearest')

        
#[FXSims,FXCcyList,dfDates] = SimulateFXRates(Path + 'FX-TimeSeries-Mod.csv',startDate,endDate,NbSims,SimLength)
#plt.plot(np.linspace(0,1,len(Sims[:,CcyList.index('JPY'),0,0])), Sims[:,CcyList.index('JPY'),0,0])

# =============================================================================
# plt.clf()
# y = []
# for i in range(0,len(Sims[:,CcyList.index('JPY'),0,0])):
#     y.append(np.average(Sims[i,CcyList.index('JPY'),:,362]))
# plt.plot(np.linspace(0,1,len(y)),y)
# =============================================================================


#FXRateIndex = dfDates.index[dfDates['DATE'] == TradeStartDate].tolist()[0] - dfDates.index[dfDates['DATE'] == startDate].tolist()[0]
a = FXfwdTrade(datetime.date(2015,6,1),datetime.date(2015,6,7),1000,'GBP',1000,'EUR')
a.GenerateMTF(datetime.date(2015,6,1),dfDates,FXCcyList,FXSims)

plt.clf()
for i in range(0,len(a.MTF[:,0])):
    plt.plot(a.EE())
    plt.plot(a.PFE(98))
    plt.plot(a.PFE(90))
    plt.plot(a.PFE(10))
    plt.plot(a.PFE(2))
plt.show()

plt.clf()
MTMVector = []
for BatchDate in daterange(a.TradeStartDate, a.maturityDate):
    a.GenerateMTF(BatchDate,dfDates,FXCcyList,FXSims)
    MTMVector.append(a.MTM())
plt.plot(MTMVector)

