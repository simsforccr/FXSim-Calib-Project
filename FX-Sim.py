# -*- coding: utf-8 -*-
"""
Soluce pour FX calib, sim and pricing

M.Jawad - 01/01/2018

"""
import datetime
import pandas as pd
import numpy as np
import scipy as sp

NbSims = 1000

class FXFwdTrade:
    startDate = None
    maturityDate = None
    
    RecLegNotional = None
    RecLegCcy = None
    RecLegCFDate = None
    
    PayLegNotional = None
    PayLegCcy = None
    PayLegCFDate = None   

def FXSim(fxspot, vol, corr, M, n, time):
    """
    FX simulation for n days with x vol
    """
    sim = np.zeros((NbSims,n))
    
    for j in range(0,n):
        if j == 0:
            sim[i,j] = fxspot
        else:
            #Sch_Rdm = np.dot(sp.linalg.cholesky(df_Corr),rdm[:,:,j-dateRange[0]])
            for i in range(0,NbSims):
                sim[i,j] = sim[0,j]*time[j]*vol*Sch_Rdm[i,j]

    return sim

def SimulateFXRates():
    """
    main function
    """
    
    Path = "C:\\Users\\Malek\\Google Drive\\FXSim-Calib-Project\\"
    startdate = datetime.date(2014,1,2)
    enddate = datetime.date(2014,12,31)
    OneYear = 252
    SimArray = []
    time_array = np.linspace(0,1,OneYear)
    
    dateparse = lambda x: pd.datetime.strptime(x, '%d/%m/%Y')
    df = pd.read_csv(Path + 'FX-TimeSeries-Mod.csv', parse_dates=['DATE'], date_parser=dateparse)
    
    # random iid standard normally distribution
    rdm = np.random.normal(0, 1, size=(6,NbSims,OneYear))

    df_LogR = np.log(df.loc[:,['AUD', 'CAD', 'EUR', 'JPY', 'CHF', 'USD']]) - np.log(df.loc[:,['AUD', 'CAD', 'EUR', 'JPY', 'CHF', 'USD']].shift(1))
    df_Vol = df_LogR.rolling(OneYear, OneYear).std()*sp.sqrt(OneYear)    
    df_Vol = pd.concat([df.loc[:,['DATE']],df_Vol], axis=1)

    dateRange = [df_Vol.index[df_Vol['DATE'] == startdate].tolist()[0], df_Vol.index[df_Vol['DATE'] == enddate].tolist()[0]]
    
    for ccy in ['AUD', 'CAD', 'EUR', 'JPY', 'CHF', 'USD']:
        ListOfSim = []
        for n in range(dateRange[0],dateRange[1]+1):
            df_Corr = df_LogR.loc[(n-OneYear):n,['AUD', 'CAD', 'EUR', 'JPY', 'CHF', 'USD']].corr(method='pearson')
            ListOfSim.append(FXSim(df.loc[n, ccy], df_Vol.loc[n,ccy], df_Corr.loc[n, ccy], rdm, OneYear, time_array))
            
        SimArray.append(ListOfSim)
    
    return SimArray

