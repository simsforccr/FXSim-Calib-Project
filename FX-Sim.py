# -*- coding: utf-8 -*-
"""
Soluce pour FX calib, sim and pricing

M.Jawad - 01/01/2018

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
NbSims = 1000
SimLength = 365
    
class FXFwdTrade:
    startDate = None
    maturityDate = None
    
    RecLegNotional = None
    RecLegCcy = None
    RecLegCFDate = None
    
    PayLegNotional = None
    PayLegCcy = None
    PayLegCFDate = None   

def FXSim(fxspot, vol, rdm, time):
    """
    FX simulation for n days with x vol using the random numbers provided over the "time" horizon
    """
    sim = np.zeros((NbSims,len(time)))
    
    for i in range(0,NbSims):
        for j in range(0,len(time)):
            sim[i,j] = fxspot+time[j]*vol*rdm[i,j]
            
    return sim

def SimulateFXRates(FilePath, SD, ED, NbSimulations, SimHorizon):
    """
    main function to simulate FX rates for each day
    input is csv file path, startdate, enddate, Number of simulations and simHorizon (nb of days)
    
    output is a 4 dimension numpy array shaped by n days, currencyIndex, Nb simulations,  simHorizon (nb of days)
        and the currency list as reference
    """
    # load historical timeseries
    dateparse = lambda x: pd.datetime.strptime(x, '%d/%m/%Y')
    df = pd.read_csv(FilePath, parse_dates=['DATE'], date_parser=dateparse)

    # setup sim variables
    SimCalibration = (ED - SD).days
    time_array = np.linspace(0,1,SimCalibration)
    CcyList = list(df)[1:]
    SimArray = np.zeros((SimCalibration,len(CcyList),NbSims,len(time_array)))
    
    df_LogR = np.log(df.loc[:,CcyList]) - np.log(df.loc[:,CcyList].shift(1))
    df_Vol = df_LogR.rolling(SimCalibration, SimCalibration).std()*sp.sqrt(SimCalibration)    
    df_Vol = pd.concat([df.loc[:,['DATE']],df_Vol], axis=1)

    dateRange = [df_Vol.index[df_Vol['DATE'] == SD].tolist()[0], df_Vol.index[df_Vol['DATE'] == ED].tolist()[0]]
    
    # generate daily FX simulations over 1 year
    for n in range(dateRange[0],dateRange[1]+1):
        df_Corr = df_LogR.loc[(n-SimCalibration):n,CcyList].corr(method='pearson')
        CorrRdm = np.zeros((6,NbSims,SimHorizon))
        # generate for every n day correlated random numbers til the SimHorizon
        for j in range(0, SimHorizon):
            CorrRdm[:,:,j] = np.dot(linalg.cholesky(df_Corr) , np.random.normal(0, 1, size=(6,NbSims)))
            
        # generate FX simulations for every n day using the Correlated Random numbers generate before
        
        for ccy in CcyList:
            ccyI = CcyList.index(ccy)
            SimArray[n-dateRange[0],ccyI,:,:] = FXSim(df.loc[n, ccy], df_Vol.loc[n,ccy], CorrRdm[ccyI,:,:], time_array)           

    return [SimArray, CcyList]

[Sims,CcyList] = SimulateFXRates(Path + 'FX-TimeSeries-Mod.csv',startDate,endDate,NbSims,SimLength)