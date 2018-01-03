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

def FXSim(fxspot, vol, rdm, Year, time):
    """
    FX simulation for n days with x vol
    """
    sim = np.zeros((NbSims,Year))
    
    for j in range(0,Year):
        for i in range(0,NbSims):
            if j == 0:
                sim[i,j] = fxspot
            else:
                sim[i,j] = sim[i,0]+time[j]*vol*rdm[i,j]
    return sim

def SimulateFXRates():
    """
    main function to simulate FX rates for each 
    """
    Path = "C:\\Users\\Malek\\Documents\\Python Projects\\FXSim-Calib-Project\\"
    #Path = "C:\\Users\\Malek\\Google Drive\\FXSim-Calib-Project\\"
    startdate = datetime.date(2014,1,2)
    enddate = datetime.date(2014,12,31)
    OneYear = 252
    SimArray = []
    time_array = np.linspace(0,1,OneYear)
    
    dateparse = lambda x: pd.datetime.strptime(x, '%d/%m/%Y')
    df = pd.read_csv(Path + 'FX-TimeSeries-Mod.csv', parse_dates=['DATE'], date_parser=dateparse)
    
    CcyList = ['AUD', 'CAD', 'EUR', 'JPY', 'CHF', 'USD']
    # random iid standard normally distribution
    rdm = np.random.normal(0, 1, size=(6,NbSims,OneYear))
    CorrRdm = np.zeros((6,NbSims,OneYear))

    df_LogR = np.log(df.loc[:,CcyList]) - np.log(df.loc[:,CcyList].shift(1))
    df_Vol = df_LogR.rolling(OneYear, OneYear).std()*sp.sqrt(OneYear)    
    df_Vol = pd.concat([df.loc[:,['DATE']],df_Vol], axis=1)

    dateRange = [df_Vol.index[df_Vol['DATE'] == startdate].tolist()[0], df_Vol.index[df_Vol['DATE'] == enddate].tolist()[0]]
    
    
    ListOfSim = []
    for n in range(dateRange[0],dateRange[1]):
        
        df_Corr = df_LogR.loc[(n-OneYear):n,CcyList].corr(method='pearson')

        rdm = np.random.normal(0, 1, size=(6,NbSims,OneYear))
        CorrRdm[:,:,n-dateRange[0]] = np.dot(linalg.cholesky(df_Corr) , rdm[:,:,n-dateRange[0]])
        
        for ccy in CcyList:
            ccyI = CcyList.index(ccy)
            ListOfSim.append(FXSim(df.loc[n, ccy], df_Vol.loc[n,ccy], CorrRdm[ccyI,:,:], OneYear, time_array))
        
    SimArray.append(ListOfSim)
    
    return SimArray